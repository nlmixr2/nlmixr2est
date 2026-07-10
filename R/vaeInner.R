# vaeInner.R -- drive the FOCEi inner likelihood directly for the VAE. The inner
# problem is set up ONCE (foceiSetup_ via vaeInnerSetup_), then the per-subject
# (and per-mixture-component) objective/gradient are evaluated through
# likInner0/lpInner in a parallel OpenMP loop (vaeInnerLik) -- reusing the inner
# logic (mixtures, multiple endpoints, error structures, log-likelihood,
# censoring) without the nlmixr2 R interface.

#' A foceiControl carrying the VAE's chosen inner likelihood + solving options.
#' focei -> interaction=1; foce -> interaction=0; laplace -> the Laplace method.
#' @noRd
.vaeInnerFoceiControl <- function(control) {
  .interaction <- if (identical(control$likelihood, "foce")) 0L else 1L
  foceiControl(rxControl = control$rxControl, maxOuterIterations = 0L,
               maxInnerIterations = 0L, covMethod = "", interaction = .interaction,
               sumProd = control$sumProd, optExpression = control$optExpression,
               literalFix = control$literalFix, literalFixRes = control$literalFixRes,
               addProp = control$addProp, calcTables = FALSE, compress = FALSE,
               eventSens = control$eventSens, indTolRelax = control$indTolRelax,
               maxOdeRecalc = control$maxOdeRecalc, odeRecalcFactor = control$odeRecalcFactor,
               stickyRecalcN = control$stickyRecalcN, print = 0L)
}

#' Set up the FOCEi inner problem for `ui` at its current ini() estimates.
#' Builds the focei opt env (which enriches the control with the model-derived
#' neta/ntheta/foceiMuGroup*/... values), augments the fit-flow-derived fields
#' (needOptimHess from the endpoint distribution, est, nF, printTop, AGQ off),
#' preprocesses the data, and calls the C++ vaeInnerSetup_ (foceiSetup_ +
#' updateTheta). Returns the setup env (keep it alive until vaeInnerFree()).
#' @noRd
.vaeInnerSetup <- function(ui, data, etaMat, control, est = "focei") {
  .ui <- rxode2::rxUiDecompress(ui)
  .fc <- .vaeInnerFoceiControl(control)
  .fc$est <- est
  .ui$control <- .fc
  .env <- .ui$foceiOptEnv
  .env$ui <- .ui
  .env$est <- est
  .env$table <- NULL
  .foceiPreProcessData(data, .env, .ui, .fc$rxControl)
  ## fit-flow-derived control fields
  .env$control$est <- est
  .env$control$printTop <- FALSE
  if (is.null(.env$control$nF)) .env$control$nF <- 0L
  .env$control$needOptimHess <- isTRUE(any(.ui$predDfFocei$distribution != "norm"))
  ## AGQ off
  .env$aqn <- 0L; .env$qx <- double(0); .env$qw <- double(0); .env$qfirst <- FALSE
  .env$nAGQ <- 0L; .env$aqLow <- -Inf; .env$aqHi <- Inf; .env$nEstOmega <- 0L
  .env$etaMat <- etaMat
  vaeInnerSetup_(.env)
  .env
}

#' Free the inner-problem state set up by .vaeInnerSetup.
#' @noRd
.vaeInnerFree <- function() invisible(vaeInnerFree_())

#' Evaluate the inner objective (and optionally the eta-gradient) at `etaMat`
#' (rows = ids: nSub, or nSub*nMix for mixtures) through the parallel C++ driver.
#' @noRd
.vaeInnerEval <- function(etaMat, control, grad = FALSE) {
  .cores <- tryCatch({
    .c <- control$rxControl$cores
    if (is.null(.c) || is.na(.c) || .c < 1L) as.integer(rxode2::getRxThreads()) else as.integer(.c)
  }, error = function(e) 1L)
  vaeInnerLik(as.matrix(etaMat), .cores, isTRUE(grad))
}
