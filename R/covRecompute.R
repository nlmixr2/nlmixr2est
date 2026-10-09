# Decoupled post-fit covariance recompute for the SAEM Louis SA-FIM ("sa") and
# the importance-sampling Monte-Carlo observed information ("imp").  Both are
# normally computed inside their own C++ kernel (src/saem.cpp, src/imp.cpp) and
# were unavailable to any other estimation method.  These helpers re-derive them
# at the converged estimates of ANY completed fit by running the native engine
# with zero estimation iterations, mirroring the pin-and-refit pattern in
# .foceiRecomputeMuCov() (R/cov.R).

#' Private decompressed copy of a fit's ui, so a nested re-fit cannot modify it
#' @param fit completed nlmixr2 fit (object or its env)
#' @return the copy, or `NULL` when it cannot be made
#' @noRd
.fitUiCopy <- function(fit) {
  tryCatch(rxode2::rxUiDecompress(unserialize(serialize(fit$ui, NULL))), error = function(e) NULL)
}

#' Pin a completed fit's converged thetas into `ui$iniDf$est`
#' @param ui rxode2 ui to pin
#' @param fit completed nlmixr2 fit
#' @return the pinned `ui`
#' @noRd
.uiPinTheta <- function(ui, fit) {
  .th <- tryCatch(fit$theta, error = function(e) NULL)
  if (!is.null(.th)) {
    .w <- match(names(.th), ui$iniDf$name)
    .ok <- !is.na(.w)
    ui$iniDf$est[.w[.ok]] <- as.numeric(.th)[.ok]
  }
  ui
}

#' A completed fit's etas as an etaMat (`$etaMat`), `NULL` when it has none
#' @param fit completed nlmixr2 fit
#' @return numeric matrix or `NULL`
#' @noRd
.fitEtaMat <- function(fit) {
  .m <- tryCatch(fit$etaMat, error = function(e) NULL)
  if (is.null(.m) || ncol(.m) == 0L) NULL else .m
}

#' Build the pinned-UI + data + etaMat needed to recompute a covariance at a
#' completed fit's converged estimates.
#'
#' The completed fit's `ui` already carries the converged theta AND omega in
#' `iniDf$est` (installed by `.nlmixr2FitUpdateParams()`); the theta is re-pinned
#' defensively.
#' @param fit completed nlmixr2 fit
#' @return list(ui, data, etaMat) or NULL on failure
#' @noRd
.covPinnedRefitArgs <- function(fit) {
  .ui <- .fitUiCopy(fit)
  if (is.null(.ui)) {
    return(NULL)
  }
  list(ui = .uiPinTheta(.ui, fit), data = getData(fit), etaMat = .fitEtaMat(fit))
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
  if (is.null(.a)) {
    return(NULL)
  }
  if (useEtaMat && !is.null(.a$etaMat)) {
    control$etaMat <- .a$etaMat
  }
  # The engines hold the parameters at the converged estimates, but they still
  # evaluate a hand-written likelihood away from them (at sampled etas, and SAEM at
  # the Omega and residual of its first iteration), where rxode2's default safeLog
  # hands back a large finite reward instead of a rejection (nlmixr2/nlmixr2est#850).
  # Ask for the log-domain mode; a no-op while every value stays valid.
  control$rxControl <- .npSafeLogDomain(control$rxControl, .a$ui)
  # the nested re-fit resets mu-referencing global state; save + restore
  .savedMuRef <- .muRefTrans$cur
  on.exit(.muRefTrans$cur <- .savedMuRef, add = TRUE)
  # fixed-parameter re-fit of an already-accepted model: bypass the prior gate
  # (#938) -- .covRecomputeFo forces est="focei", which declares no prior
  # support, and the try() below would otherwise silently return NULL for a
  # prior-carrying fit
  .fit2 <- try(
    suppressMessages(suppressWarnings(
      .nlmixr2PriorGateBypass(
        nlmixr2(.a$ui, data = .a$data, est = est, control = control)
      )
    )),
    silent = TRUE
  )
  if (inherits(.fit2, "try-error")) {
    return(NULL)
  }
  .cov <- tryCatch(.fit2$cov, error = function(e) NULL)
  if (is.null(.cov) || !is.matrix(.cov)) {
    return(NULL)
  }
  # .fit2 is a full re-fit, so its $cov has already been through
  # .mixInstallProbScaleCov(); say so, or .covInstallResult() rotates it twice.
  list(cov = .cov, covMethod = .fit2$covMethod, mixRotated = TRUE)
}

#' Recompute the SAEM Louis SA-FIM ("sa") at any fit's converged estimates.
#'
#' Runs a short SAEM started at the pinned (converged) theta/omega, with every
#' population parameter held there (`saemHoldPar`, `.saemHoldCfg()`): the
#' `nBurn`/`nEm` warm-up only equilibrates the MCMC chains, and the dedicated
#' `nSaCov` phase accumulates the Louis observed information at the fit's own
#' estimates.
#' @param fit completed nlmixr2 fit
#' @param control `saControl()` options, or `NULL` for the defaults
#' @return list(cov, covMethod, extras) or NULL
#' @noRd
.covRecomputeSa <- function(fit, control = NULL) {
  if (is.null(control)) {
    control <- saControl()
  }
  .state <- if (control$warmStart) .saemChainState(fit)
  # SAEM derives its own etaMat from the MCMC; no external eta seed
  .covRecomputeNative(fit, "saem", .covEngineControl("sa", control, .state), useEtaMat = FALSE)
}

#' The state a SAEM fit's covariance phase continues from
#' @param fit nlmixr2 fit
#' @return `NULL` when the fit kept no chains, otherwise a list: `phiM`, the
#'   `(N * nmc) x nphi` chain state at the last estimation iteration;
#'   `sigma2`, the per-endpoint sigma2 the fit's Louis residual score last
#'   read; and `mpostPhi`, the `N x nphi` posterior means (each `NULL` when the
#'   fit does not have it)
#' @noRd
.saemChainState <- function(fit) {
  .phiM <- tryCatch(fit$phiM, error = function(e) NULL)
  if (!is.array(.phiM) || length(dim(.phiM)) != 4L || any(dim(.phiM) == 0L)) {
    return(NULL)
  }
  .d <- dim(.phiM)
  .last <- .phiM[,, .d[3], , drop = FALSE]
  dim(.last) <- c(.d[1] * .d[2], .d[4])
  if (anyNA(.last)) {
    return(NULL)
  }
  .saem <- tryCatch(fit$saem, error = function(e) NULL)
  .sigma2 <- tryCatch(as.numeric(.saem$res_info$sigma2), error = function(e) NULL)
  if (length(.sigma2) == 0L || !all(is.finite(.sigma2)) || any(.sigma2 <= 0)) {
    .sigma2 <- NULL
  }
  .mpost <- tryCatch(.saem$mpost_phi, error = function(e) NULL)
  if (!is.matrix(.mpost) || !identical(dim(.mpost), .d[c(1L, 4L)]) || !all(is.finite(.mpost))) {
    .mpost <- NULL
  }
  list(phiM = .last, sigma2 = .sigma2, mpostPhi = .mpost)
}

#' Recompute the importance-sampling Monte-Carlo covariance ("imp") at any fit's
#' converged estimates.
#'
#' Runs the imp kernel at the pinned converged estimates with `impFrozen`:
#' `nIter` E-steps (`mapIter=0`) and no M-step, so the parameters stay where
#' they are and the MAP pass + `impComputeCov` evaluate the Monte-Carlo
#' observed information at the fit's own estimates.
#' @param fit completed nlmixr2 fit
#' @param control `impCovControl()` options, or `NULL` for the defaults
#' @return list(cov, covMethod, extras) or NULL
#' @noRd
.covRecomputeImp <- function(fit, control = NULL) {
  .covRecomputeNative(fit, "imp", .covEngineControl("imp", control), useEtaMat = TRUE)
}

#' Engine control for a decoupled covariance recompute
#' @param method "sa" or "imp"
#' @param control `saControl()`/`impCovControl()` options, or `NULL` for the
#'   defaults
#' @param state for "sa", a SAEM fit's chain state (`.saemChainState()`) to
#'   continue with no warm-up iterations, or `NULL` to start the chains around
#'   the estimates and run `nBurn`/`nEm` warm-up iterations
#' @return `saemControl()` or `impmapControl()` object
#' @noRd
.covEngineControl <- function(method, control = NULL, state = NULL) {
  if (identical(method, "sa")) {
    if (is.null(control)) {
      control <- saControl()
    }
    # a SAEM fit's own iterations are the warm-up of its chains
    .warm <- !is.null(state)
    return(saemControl(
      nBurn = if (.warm) 0L else control$nBurn,
      nEm = if (.warm) 0L else control$nEm,
      nSaCov = control$nSaCov,
      seed = control$seed,
      covMethod = "sa",
      calcTables = FALSE,
      saemHoldPar = TRUE,
      saemWarmState = state
    ))
  }
  if (is.null(control)) {
    control <- impCovControl()
  }
  # impmap's default SIR sample (at least 25) cannot exceed a small isample
  .sir <- min(max(25L, as.integer(ceiling(max(control$isample) / 10))), min(control$isample))
  impmapControl(
    nIter = control$nIter,
    mapIter = 0L,
    isample = control$isample,
    impSeed = control$impSeed,
    sirSample = .sir,
    covMethod = "imp",
    calcTables = FALSE,
    impFrozen = TRUE
  )
}

#' Dispatcher: recompute a decoupled covariance ("sa"/"imp") on a completed fit.
#' @param fit completed nlmixr2 fit
#' @param method "sa" or "imp"
#' @param control covariance control, or `NULL` for the defaults
#' @return list(cov, covMethod, extras) or NULL
#' @noRd
.covRecompute <- function(fit, method, control = NULL) {
  if (identical(method, "sa")) {
    return(.covRecomputeSa(fit, control = control))
  }
  if (identical(method, "imp")) {
    return(.covRecomputeImp(fit, control = control))
  }
  NULL
}

#' Install a recompute result (list(cov, covMethod, mixRotated)) onto a fit env.
#'
#' Installs through `.covInstall()`: a matrix `.covGuard()` rejects is NOT
#' installed (the existing covariance is kept, never silently downgraded), the
#' prior covariance stays recoverable via `covList`/`setCov()`, and SE/%RSE/CI
#' are refreshed on the fit's OWN parameter table (its point estimates are kept).
#' @param env fit environment
#' @param r recompute result from `.covRecompute()` (or NULL)
#' @param warn warn when nothing usable was installed (see `.covInstall()`)
#' @param what requested covariance-method name, for the warnings
#' @return invisibly TRUE if a new covariance was installed
#' @noRd
.covInstallResult <- function(env, r, warn = FALSE, what = r$covMethod) {
  .cov <- r$cov
  # A covariance computed directly (analytic) is on the mlogit estimation scale
  # and needs the mixture block rotated onto the probability scale; one that came
  # back from a re-fit (sa/imp, via .covRecompute) was already rotated there.
  if (!isTRUE(r$mixRotated)) {
    .cov <- .covToReportedScale(env, .cov)
  }
  # An engine whose covariance carries no mixture rows at all (saem) needs the
  # (7.51) block appended again -- the recomputed matrix REPLACED the one that
  # had it, so without this a setCov() drops the proportion back to SE = NA.
  .covInstall(env, .cov, r$covMethod, what = what, warn = warn, mixAppend = TRUE)
}

#' Read the deferred foreign-covariance request ("sa"/"imp") stashed on a fit's
#' control by the control resolver.
#' @param fit completed nlmixr2 fit (object or env)
#' @return "sa"/"imp" or NA_character_
#' @noRd
.covGetDeferred <- function(fit) {
  # fit$control is the uniform per-method control accessor (nmObjGetControl);
  # families that finalize through .foceiFamilyReturn (vae/vi/impmap/np) carry
  # the deferred request on the internal foceiControl instead.
  for (.acc in c("control", "foceiControl")) {
    .ctl <- tryCatch(do.call("$", list(fit, .acc)), error = function(e) NULL)
    .d <- if (is.list(.ctl)) .ctl$covMethodDeferred else NULL
    if (!is.null(.d) && length(.d) == 1L && !is.na(.d) && nzchar(.d)) return(.d)
  }
  NA_character_
}
