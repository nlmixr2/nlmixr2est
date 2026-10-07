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
  # This re-fit is pinned at the converged estimates but still takes a frozen EM / SA
  # step, so it CAN move a theta.  For a hand-written likelihood that means it can step a
  # scale out of its domain, where rxode2's default safeLog hands back a large finite
  # reward instead of a rejection -- and the covariance would then be formed around a
  # point the likelihood cannot evaluate (nlmixr2/nlmixr2est#850).  Ask for the
  # log-domain mode here too; a no-op while every parameter stays valid.
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
#' Runs a short SAEM at the pinned (converged) theta/omega -- a modest warm-up
#' (`nBurn`/`nEm`) equilibrates the MCMC chains and the stochastic-approximation
#' running sums (a cold `nBurn=0, nEm=0` start leaves those uninitialized and
#' produces a non-finite FIM), then the dedicated `nSaCov` phase accumulates the
#' Louis observed-information at the (essentially unchanged) converged point.
#' @param fit completed nlmixr2 fit
#' @param control `saControl()` options, or `NULL` for the defaults
#' @return list(cov, covMethod, extras) or NULL
#' @noRd
.covRecomputeSa <- function(fit, control = NULL) {
  # SAEM derives its own etaMat from the MCMC; no external eta seed
  .covRecomputeNative(fit, "saem", .covEngineControl("sa", control), useEtaMat = FALSE)
}

#' Recompute the importance-sampling Monte-Carlo covariance ("imp") at any fit's
#' converged estimates.
#'
#' Runs the impmap kernel (already `maxOuterIterations=0`) with a single frozen
#' EM step (`nIter=1, mapIter=0`) at the
#' pinned converged estimates, so the MAP pass + `impComputeCov` evaluate the
#' Monte-Carlo observed information essentially at the converged point.
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
#' @return `saemControl()` or `impmapControl()` object
#' @noRd
.covEngineControl <- function(method, control = NULL) {
  if (identical(method, "sa")) {
    if (is.null(control)) {
      control <- saControl()
    }
    return(saemControl(
      nBurn = control$nBurn,
      nEm = control$nEm,
      nSaCov = control$nSaCov,
      seed = control$seed,
      covMethod = "sa",
      calcTables = FALSE
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
    calcTables = FALSE
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
