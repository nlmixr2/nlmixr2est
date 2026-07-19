# Inner-problem harness for the nonparametric engines.  Sets up the FOCEi inner
# problem once (via vaeInnerSetup_) so the conditional-likelihood primitive
# (npEvalCondLik) and the Psi builder (npBuildPsi) can be evaluated per support
# point, reusing the ODE solve, residual-error models, transform-both-sides and
# censoring unchanged.  Mirrors .fsaemInnerSetup / .adviInnerSetup.

#' A foceiControl carrying the npag/npb inner likelihood + solving options.
#' @noRd
.npInnerFoceiControl <- function(control) {
  foceiControl(rxControl = control$rxControl, maxOuterIterations = 0L,
               maxInnerIterations = 0L, covMethod = "", interaction = 1L,
               sumProd = control$sumProd, optExpression = control$optExpression,
               literalFix = control$literalFix, literalFixRes = control$literalFixRes,
               addProp = control$addProp, calcTables = FALSE, compress = FALSE,
               maxOdeRecalc = control$maxOdeRecalc, odeRecalcFactor = control$odeRecalcFactor,
               stickyRecalcN = control$stickyRecalcN, print = 0L)
}

#' Set up the FOCEi inner problem for the nonparametric engines.
#' @param ui rxode2 ui object (already bounded-transformed by the dispatch hook)
#' @param data estimation data
#' @param etaMat starting etas [nsub, neta] (support points are supplied later)
#' @param control an impmapControl-derived control
#' @return the setup env (keep alive until .npInnerFree())
#' @noRd
.npInnerSetup <- function(ui, data, etaMat, control) {
  .ui <- rxode2::rxUiDecompress(ui)
  .fc <- .npInnerFoceiControl(control)
  .fc$est <- "focei"
  .ui$control <- .fc
  .env <- .ui$foceiOptEnv
  .env$ui <- .ui
  .env$est <- "focei"
  .env$table <- NULL
  .foceiPreProcessData(data, .env, .ui, .fc$rxControl)
  .env$control$est <- "focei"
  .env$control$printTop <- FALSE
  if (is.null(.env$control$nF)) .env$control$nF <- 0L
  .env$control$needOptimHess <- isTRUE(any(.ui$predDfFocei$distribution != "norm"))
  .env$aqn <- 0L; .env$qx <- double(0); .env$qw <- double(0); .env$qfirst <- FALSE
  .env$nAGQ <- 0L; .env$aqLow <- -Inf; .env$aqHi <- Inf; .env$nEstOmega <- 0L
  .env$etaMat <- etaMat
  vaeInnerSetup_(.env)
  .env
}

#' Free the inner-problem state set up by .npInnerSetup.
#' @noRd
.npInnerFree <- function() invisible(vaeInnerFree_())

#' Build the Psi conditional-likelihood matrix (subjects x support points) at the
#' supplied eta support points, on the already set-up inner problem.
#' @param etaPoints matrix of support points, one per row (columns = etas)
#' @param control an impmapControl-derived control (for the thread count)
#' @return numeric matrix psi (subjects x support points)
#' @noRd
.npInnerPsi <- function(etaPoints, control) {
  .cores <- tryCatch({
    .c <- control$rxControl$cores
    if (is.null(.c) || is.na(.c) || .c < 1L) as.integer(rxode2::getRxThreads()) else as.integer(.c)
  }, error = function(e) 1L)
  npBuildPsi(as.matrix(etaPoints), .cores)
}
