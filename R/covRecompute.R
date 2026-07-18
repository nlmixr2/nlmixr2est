# Decoupled post-fit covariance recompute for the SAEM Louis SA-FIM ("sa") and
# the importance-sampling Monte-Carlo observed information ("imp").  Both are
# normally computed inside their own C++ kernel (src/saem.cpp, src/imp.cpp) and
# were unavailable to any other estimation method.  These helpers re-derive them
# at the converged estimates of ANY completed fit by running the native engine
# with zero estimation iterations, mirroring the pin-and-refit pattern in
# .foceiRecomputeMuCov() (R/cov.R).

#' Build the pinned-UI + data + etaMat needed to recompute a covariance at a
#' completed fit's converged estimates.
#'
#' The completed fit's `ui` already carries the converged theta AND omega in
#' `iniDf$est` (installed by `.nlmixr2FitUpdateParams()`); the theta is re-pinned
#' defensively.  A deep copy is returned so the nested re-fit cannot mutate this
#' fit's UI.
#' @param fit completed nlmixr2 fit
#' @return list(ui, data, etaMat) or NULL on failure
#' @noRd
.covPinnedRefitArgs <- function(fit) {
  .ui <- tryCatch(rxode2::rxUiDecompress(unserialize(serialize(fit$ui, NULL))),
                  error = function(e) NULL)
  if (is.null(.ui)) return(NULL)
  .th <- tryCatch(fit$theta, error = function(e) NULL)
  if (!is.null(.th)) {
    .w <- match(names(.th), .ui$iniDf$name)
    .ok <- !is.na(.w)
    .ui$iniDf$est[.w[.ok]] <- as.numeric(.th)[.ok]
  }
  .eta <- tryCatch(fit$eta, error = function(e) NULL)
  .etaMat <- NULL
  if (!is.null(.eta)) {
    .etaCols <- setdiff(names(.eta), "ID")
    if (length(.etaCols) > 0L) {
      .etaMat <- as.matrix(.eta[, .etaCols, drop = FALSE])
    }
  }
  list(ui = .ui, data = getData(fit), etaMat = .etaMat)
}

#' Run a native engine (saem/imp) at the pinned converged estimates and harvest
#' its covariance + already-rendered parameter table.
#' @param fit completed nlmixr2 fit
#' @param est native engine to run ("saem" or "imp")
#' @param control control object for that engine (zero-iteration + cov request)
#' @param useEtaMat whether the engine accepts an `etaMat` seed
#' @return list(cov, covMethod, extras) or NULL
#' @noRd
.covRecomputeNative <- function(fit, est, control, useEtaMat = TRUE) {
  .a <- .covPinnedRefitArgs(fit)
  if (is.null(.a)) return(NULL)
  if (useEtaMat && !is.null(.a$etaMat)) control$etaMat <- .a$etaMat
  # the nested re-fit resets mu-referencing global state; save + restore
  .savedMuRef <- .muRefTrans$cur
  on.exit(.muRefTrans$cur <- .savedMuRef, add = TRUE)
  .fit2 <- try(suppressMessages(suppressWarnings(
    nlmixr2(.a$ui, data = .a$data, est = est, control = control))), silent = TRUE)
  if (inherits(.fit2, "try-error")) return(NULL)
  .cov <- tryCatch(.fit2$cov, error = function(e) NULL)
  if (is.null(.cov) || !is.matrix(.cov)) return(NULL)
  list(cov = .cov, covMethod = .fit2$covMethod)
}

#' Recompute the SAEM Louis SA-FIM ("sa") at any fit's converged estimates.
#'
#' Runs a short SAEM at the pinned (converged) theta/omega -- a modest warm-up
#' (`nBurn`/`nEm`) equilibrates the MCMC chains and the stochastic-approximation
#' running sums (a cold `nBurn=0, nEm=0` start leaves those uninitialized and
#' produces a non-finite FIM), then the dedicated `nSaCov` phase accumulates the
#' Louis observed-information at the (essentially unchanged) converged point.
#' @param fit completed nlmixr2 fit
#' @param nBurn,nEm short warm-up iteration counts (default 100 each)
#' @return list(cov, covMethod, extras) or NULL
#' @noRd
.covRecomputeSa <- function(fit, nBurn = 100L, nEm = 100L) {
  # SAEM derives its own etaMat from the MCMC; no external eta seed
  .covRecomputeNative(fit, "saem",
                      saemControl(nBurn = as.integer(nBurn), nEm = as.integer(nEm),
                                  covMethod = "sa", calcTables = FALSE),
                      useEtaMat = FALSE)
}

#' Recompute the importance-sampling Monte-Carlo covariance ("imp") at any fit's
#' converged estimates.
#'
#' Runs the impmap kernel (already `maxOuterIterations=0`) with a single frozen
#' EM step (`nIter=1, mapIter=0`; the kernel requires >=1 iteration) at the
#' pinned converged estimates, so the MAP pass + `impComputeCov` evaluate the
#' Monte-Carlo observed information essentially at the converged point.
#' @param fit completed nlmixr2 fit
#' @param nIter frozen EM iterations (default 1; kernel segfaults at 0)
#' @return list(cov, covMethod, extras) or NULL
#' @noRd
.covRecomputeImp <- function(fit, nIter = 1L) {
  .covRecomputeNative(fit, "imp",
                      impmapControl(nIter = as.integer(nIter), mapIter = 0L,
                                    covMethod = "imp", calcTables = FALSE),
                      useEtaMat = TRUE)
}

#' Dispatcher: recompute a decoupled covariance ("sa"/"imp") on a completed fit.
#' @param fit completed nlmixr2 fit
#' @param method "sa" or "imp"
#' @return list(cov, covMethod, extras) or NULL
#' @noRd
.covRecompute <- function(fit, method) {
  if (identical(method, "sa")) return(.covRecomputeSa(fit))
  if (identical(method, "imp")) return(.covRecomputeImp(fit))
  NULL
}

#' Install a recompute result (list(cov, covMethod)) onto a fit env.
#'
#' PD-guards the incoming covariance (a non-finite / non-PD matrix is NOT
#' installed -- the existing covariance is kept, never silently downgraded),
#' keeps the prior covariance recoverable via `covList`/`setCov()`, and refreshes
#' SE/%RSE/CI on the fit's OWN parameter table from the new covariance (the base
#' fit's point estimates are preserved).  Mirrors `.saemInstallAnalyticCov()`.
#' @param env fit environment
#' @param r recompute result from `.covRecompute()` (or NULL)
#' @return invisibly TRUE if a new covariance was installed
#' @noRd
.covInstallResult <- function(env, r) {
  if (is.null(r) || is.null(r$cov) || !is.matrix(r$cov)) return(invisible(FALSE))
  .cov <- 0.5 * (r$cov + t(r$cov))                         # exact symmetry
  .ev <- suppressWarnings(eigen(.cov, symmetric = TRUE, only.values = TRUE)$values)
  if (any(!is.finite(diag(.cov))) || any(diag(.cov) <= 0) ||
        !all(is.finite(.ev)) || min(.ev) <= 0) {
    return(invisible(FALSE))                               # keep the existing cov
  }
  # keep the prior covariance recoverable via setCov()
  if (exists("cov", envir = env, inherits = FALSE) && is.matrix(env$cov)) {
    .stash <- list(env$cov)
    names(.stash) <- as.character(if (exists("covMethod", envir = env, inherits = FALSE)) {
      env$covMethod
    } else "prev")
    .cl <- if (exists("covList", envir = env, inherits = FALSE)) env$covList else NULL
    if (is.null(.cl[[names(.stash)]]) && !identical(names(.stash), r$covMethod)) {
      .cl <- c(.cl, .stash)
    }
    assign("covList", .cl, envir = env)
  }
  assign("cov", .cov, envir = env)
  assign("covMethod", r$covMethod, envir = env)
  # refresh SE/%RSE/CI on the fit's own parameter table from the new covariance
  if (exists("parFixedDf", envir = env, inherits = FALSE)) {
    .pf <- env$parFixedDf
    .se <- sqrt(diag(.cov))
    .ci <- tryCatch(as.numeric(rxode2::rxGetControl(env$ui, "ci", 0.95)),
                    error = function(e) 0.95)
    .qn <- stats::qnorm(1 - (1 - .ci) / 2)
    for (.n in rownames(.pf)) {
      if (.n %in% names(.se) && "SE" %in% names(.pf)) {
        .s <- .se[[.n]]; .e <- .pf[.n, "Estimate"]
        .pf[.n, "SE"] <- .s
        if ("%RSE" %in% names(.pf)) .pf[.n, "%RSE"] <- abs(.s / .e) * 100
        if (all(c("CI Lower", "CI Upper", "Back-transformed") %in% names(.pf)) &&
              isTRUE(all.equal(unname(.pf[.n, "Back-transformed"]), unname(.e)))) {
          .pf[.n, "CI Lower"] <- .e - .qn * .s
          .pf[.n, "CI Upper"] <- .e + .qn * .s
        }
      }
    }
    env$parFixedDf <- .pf
  }
  .nlmixr2CovConditionUpdate(env)
  invisible(TRUE)
}

#' Read the deferred foreign-covariance request ("sa"/"imp") stashed on a fit's
#' control by the control resolver.
#' @param fit completed nlmixr2 fit (object or env)
#' @return "sa"/"imp" or NA_character_
#' @noRd
.covGetDeferred <- function(fit) {
  # the family controls are reconstructed through the nmObjGet `$` accessor, not
  # stored as literal env variables -- read them via `$` (do.call keeps the name
  # a variable) rather than get() on the fit env.
  for (.cn in c("foceiControl", "saemControl", "nlmControl", "vaeControl",
                "adviControl", "nlmeControl")) {
    .ctl <- tryCatch(do.call("$", list(fit, .cn)), error = function(e) NULL)
    .d <- if (is.list(.ctl)) .ctl$covMethodDeferred else NULL
    if (!is.null(.d) && length(.d) == 1L && !is.na(.d) && nzchar(.d)) return(.d)
  }
  NA_character_
}
