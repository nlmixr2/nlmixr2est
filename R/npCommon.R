# Shared setup for the nonparametric estimation methods (npag, npb) and their
# mu-referenced sugar (mnpag/inpag, mnpb/inpb).  All route through the shared
# FOCEI family plumbing (reusing the impmap family fit) so the inner model, mu
# index maps, covariance, and tables are inherited; the C++ driver (npagOuter /
# npbOuter) is selected by op_focei.isNpag / isNpb from the est string.
#
# `muModel` is NULL for the plain methods (respect the user control) or forced to
# "lin" (OLS covariate M-step) / "irls" (reweighted) for the sugar variants.  The
# lin-vs-irls covariate M-step itself is implemented in a later milestone; here
# the sugar records the intent and selects the same driver.

#' Reject generalized (non-normal) likelihoods for the nonparametric engines.
#' The npag/npb Psi is the conditional density from the FOCEi inner problem; the
#' residual-error (gamma) handling and the -2LL only make sense for normally-
#' distributed endpoints.  Censoring (BLQ/ALQ) and transform-both-sides stay
#' "norm" and are allowed; discrete / user-`ll()` endpoints are not.
#' @noRd
.npAssertNormal <- function(ui, est) {
  .dist <- tryCatch(ui$predDfFocei$distribution, error = function(e) NULL)
  if (is.null(.dist)) {
    .dist <- tryCatch(ui$predDf$distribution, error = function(e) NULL)
  }
  if (!is.null(.dist) && length(.dist) > 0L && any(as.character(.dist) != "norm")) {
    stop("the '", est, "' estimation routine does not support generalized ",
         "(non-normal) likelihoods; only normally-distributed endpoints -- ",
         "optionally with censoring (BLQ/ALQ) or transform-both-sides -- are ",
         "supported", call. = FALSE)
  }
  invisible()
}

#' @noRd
.npEstCore <- function(env, est, muModel = NULL, ...) {
  .ui <- env$ui
  .what <- paste0(" for the estimation routine '", est, "'")
  .npAssertNormal(.ui, est)
  if (!rxode2hasLlik()) {
    rxode2::assertRxUiTransformNormal(.ui, .what, .var.name = .ui$modelName)
  }
  rxode2::assertRxUiIovNoCor(.ui, .what, .var.name = .ui$modelName)
  .foceiFamilyControl(env, ..., type = "impmapControl")
  # .foceiFamilyControl populated ui$control with the full foceiControl fields
  # (needOptimHess, etc.); read that authoritative control, add the nonparametric
  # box + knobs to it, and write it back to both ui$control (read by
  # .impmapFamilyFit) and env$control so npagOuter (src/npag.cpp) sees them in
  # e$control while the standard foceiControl fields are preserved.
  .ctl <- get("control", envir = .ui)
  if (!is.null(muModel)) {          # sugar variant: force mu-referencing on
    .ctl$muModel <- muModel
    .ctl$muRefCovAlg <- TRUE
  }
  .box <- .npEtaBox(.ui, .ctl)
  .ctl$npBoxLower <- as.numeric(.box$lower)
  .ctl$npBoxUpper <- as.numeric(.box$upper)
  .ctl$npEtaNames <- .box$names
  .ctl$npPoints <- as.integer(if (is.null(.ctl$points)) 2028L else .ctl$points)
  .ctl$npCycles <- as.integer(if (is.null(.ctl$cycles)) 100L else .ctl$cycles)
  .ctl$npGammaOptimize <-
    isTRUE(if (is.null(.ctl$gammaOptimize)) TRUE else .ctl$gammaOptimize)
  assign("control", .ctl, envir = .ui)
  env$control <- .ctl
  .control <- .ctl
  on.exit({
    if (is.environment(.ui) && exists("control", envir = .ui, inherits = FALSE)) {
      rm("control", envir = .ui)
    }
  }, add = TRUE)
  env$impmapControl <- .control
  env$est <- est
  .ui <- env$ui
  .npFamilyFit(env, .ui, ...)
}

# Fit driver for the nonparametric engines.  Turns off the FOCEI outer optimizer
# (npagOuter/npbOuter drive the cycle), builds the 0-based mu index maps that the
# finalization (impMuInterceptStep) reuses, and -- unlike .impmapFamilyFit --
# does NOT build or wire the theta-sensitivity model (npag/npb evaluate the
# conditional likelihood directly; the sens-augmented inner model breaks the
# per-observation offsets that npBuildPsiCore relies on).
#' @noRd
.npFamilyFit <- function(env, ui, ...) {
  .control <- ui$control
  .control$maxOuterIterations <- 0L
  .control$covMethod <- 0L
  .env <- ui$foceiOptEnv     # builds foceiMuGroupTheta (covariate mu-groups)
  .iniDf <- ui$iniDf
  .th <- .iniDf[!is.na(.iniDf$ntheta), ]
  .thNames <- .th[order(.th$ntheta), "name"]
  .etaRows <- .iniDf[!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2, ]
  .etaNames <- .etaRows[order(.etaRows$neta1), "name"]
  .mr <- ui$muRefDataFrame
  .muThetaIdx <- as.integer(match(.mr$theta, .thNames) - 1L)
  .muEtaIdx <- as.integer(match(.mr$eta, .etaNames) - 1L)
  .covGroupTheta <- rxode2::rxGetControl(ui, "foceiMuGroupTheta", integer(0))
  .keep <- !is.na(.muThetaIdx) & !is.na(.muEtaIdx) & !(.muThetaIdx %in% .covGroupTheta)
  .control$impMuThetaIdx <- .muThetaIdx[.keep]
  .control$impMuEtaIdx <- .muEtaIdx[.keep]
  .control$impThetaSensIdx <- integer(0)   # no sensitivity model for npag/npb
  .etaOrd <- .etaRows[order(.etaRows$neta1), ]
  .control$impOmegaFixedEta <- as.integer(which(isTRUE(.etaOrd$fix) | .etaOrd$fix) - 1L)
  assign("control", .control, envir = ui)
  .est <- if (exists("est", envir = env)) get("est", envir = env) else "npag"
  .foceiFamilyReturn(env, ui, ..., est = .est)
}

# mu-attribute for the plain methods: gated on the control (mu only when the user
# asked for it), matching the imp / mfocei families.
.npMuAttr <- function(control) {
  isTRUE(!identical(control$muModel, "none")) && isTRUE(control$muRefCovAlg)
}

# Validate a control for a nonparametric engine.  The impmap validator rebuilds
# the control via do.call(impmapControl, .), which rejects the npag-only fields
# (points/cycles/gammaOptimize/est), so strip them first, then re-attach.
#' @noRd
.npValidCtl <- function(control, est) {
  .in <- control[[1]]
  .pts <- if (is.null(.in$points)) 2028L else .in$points
  .cyc <- if (is.null(.in$cycles)) 100L else .in$cycles
  .go  <- if (is.null(.in$gammaOptimize)) TRUE else .in$gammaOptimize
  if (is.list(.in)) {
    .in$points <- NULL; .in$cycles <- NULL; .in$gammaOptimize <- NULL; .in$est <- NULL
  }
  .ctl <- getValidNlmixrCtl.impmap(list(.in))
  .ctl$est <- est
  .ctl$points <- as.integer(.pts)
  .ctl$cycles <- as.integer(.cyc)
  .ctl$gammaOptimize <- isTRUE(.go)
  .ctl
}

# Control validators / accessors for the nonparametric methods and their sugar.
# All delegate to the impmap family (the shared FOCEI-family control) and stamp
# the est string so the C++ driver dispatch is preserved.

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mnpag <- function(control) {
  .ctl <- .npValidCtl(control, "mnpag"); .ctl$muModel <- "lin"; .ctl
}
#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.inpag <- function(control) {
  .ctl <- .npValidCtl(control, "inpag"); .ctl$muModel <- "irls"; .ctl
}
#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mnpb <- function(control) {
  .ctl <- .npValidCtl(control, "mnpb"); .ctl$muModel <- "lin"; .ctl
}
#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.inpb <- function(control) {
  .ctl <- .npValidCtl(control, "inpb"); .ctl$muModel <- "irls"; .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.mnpag <- function(x, ...) nmObjGetControl.impmap(x, ...)
#' @rdname nmObjGetControl
#' @export
nmObjGetControl.inpag <- function(x, ...) nmObjGetControl.impmap(x, ...)
#' @rdname nmObjGetControl
#' @export
nmObjGetControl.mnpb <- function(x, ...) nmObjGetControl.impmap(x, ...)
#' @rdname nmObjGetControl
#' @export
nmObjGetControl.inpb <- function(x, ...) nmObjGetControl.impmap(x, ...)

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.mnpag <- function(x, ...) nmObjGetFoceiControl.impmap(x, ...)
#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.inpag <- function(x, ...) nmObjGetFoceiControl.impmap(x, ...)
#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.mnpb <- function(x, ...) nmObjGetFoceiControl.impmap(x, ...)
#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.inpb <- function(x, ...) nmObjGetFoceiControl.impmap(x, ...)
