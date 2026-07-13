# advi.R -- orchestration for est="advi" (automatic differentiation variational
# inference, Kucukelbir et al. 2017).  Sets up the FOCEi inner problem (reused
# for the per-subject log-joint and eta-gradient) plus, when non-mu structural
# thetas are present, the impmap theta-sensitivity model (reused for the outer
# population gradient), then drives the ADVI optimization loop in C++.

#' A foceiControl carrying ADVI's chosen inner likelihood + solving options.
#' Mirrors .vaeInnerFoceiControl: focei -> interaction=1; foce/focep ->
#' interaction=0 (focep = FOCE+, R at the live conditional eta); laplace -> the
#' Laplace method.
#' @noRd
.adviInnerFoceiControl <- function(control) {
  .lik <- control$likelihood
  .interaction <- if (.lik %in% c("foce", "focep")) 0L else 1L
  .foce <- if (identical(.lik, "focep")) "foce+" else "nonmem"
  foceiControl(rxControl = control$rxControl, maxOuterIterations = 0L,
               maxInnerIterations = 0L, covMethod = "", interaction = .interaction,
               foce = .foce,
               sumProd = control$sumProd, optExpression = control$optExpression,
               literalFix = control$literalFix, literalFixRes = control$literalFixRes,
               addProp = control$addProp, calcTables = FALSE, compress = FALSE,
               eventSens = control$eventSens, indTolRelax = control$indTolRelax,
               maxOdeRecalc = control$maxOdeRecalc, odeRecalcFactor = control$odeRecalcFactor,
               stickyRecalcN = control$stickyRecalcN, print = 0L)
}

#' Set up the FOCEi inner problem (reused for the per-subject log-joint and
#' eta-gradient) plus, when non-mu structural/sigma thetas are present, the
#' impmap theta-sensitivity model (reused for the outer population gradient).
#' Modeled on .vaeInnerSetup, adding the 0-based `impThetaSensIdx` so foceiSetup_
#' wires the sensitivity output offsets into op_focei.
#' @param ui rxode2 ui object (already bounded-transformed by the dispatch hook)
#' @param data estimation data
#' @param etaMat starting etas [nsub, neta]
#' @param control adviControl
#' @return the setup env (keep alive until .adviInnerFree())
#' @noRd
.adviInnerSetup <- function(ui, data, etaMat, control) {
  .ui <- rxode2::rxUiDecompress(ui)
  .fc <- .adviInnerFoceiControl(control)
  .fc$est <- "advi"
  ## 0-based non-mu theta indices with d(f)/d(theta) & d(V)/d(theta) outputs; the
  ## theta-sensitivity model (built in .foceiOptEnvLik for est="advi") supplies
  ## the columns and foceiSetup_ records their lhs offsets in op_focei.
  .fc$impThetaSensIdx <- as.integer(.impmapEstTheta(.ui)$all - 1L)
  .ui$control <- .fc
  .env <- .ui$foceiOptEnv
  .env$ui <- .ui
  .env$est <- "advi"
  .env$table <- NULL
  .foceiPreProcessData(data, .env, .ui, .fc$rxControl)
  .env$control$est <- "advi"
  ## foceiSetup_ reads impThetaSensIdx from e$control (foceiO); make sure it is
  ## present there (not only on the pre-build .fc) so op_focei wires the offsets.
  .env$control$impThetaSensIdx <- as.integer(.impmapEstTheta(.ui)$all - 1L)
  .env$control$printTop <- FALSE
  if (is.null(.env$control$nF)) .env$control$nF <- 0L
  .env$control$needOptimHess <- isTRUE(any(.ui$predDfFocei$distribution != "norm"))
  .env$aqn <- 0L; .env$qx <- double(0); .env$qw <- double(0); .env$qfirst <- FALSE
  .env$nAGQ <- 0L; .env$aqLow <- -Inf; .env$aqHi <- Inf; .env$nEstOmega <- 0L
  .env$etaMat <- etaMat
  vaeInnerSetup_(.env)
  .env
}

#' Free the inner-problem state set up by .adviInnerSetup.
#' @noRd
.adviInnerFree <- function() invisible(vaeInnerFree_())

#' Evaluate the inner objective (and optionally the eta-gradient) at `etaMat`
#' (rows = ids) through the parallel C++ driver (reused verbatim from vae).
#' @noRd
.adviInnerEval <- function(etaMat, control, grad = FALSE, preds = FALSE) {
  .cores <- tryCatch({
    .c <- control$rxControl$cores
    if (is.null(.c) || is.na(.c) || .c < 1L) as.integer(rxode2::getRxThreads()) else as.integer(.c)
  }, error = function(e) 1L)
  vaeInnerLik(as.matrix(etaMat), .cores, isTRUE(grad), isTRUE(preds))
}

#' Fit an ADVI model: set up the inner/outer problems and run the C++ loop.
#' @param env estimation environment (holds ui, data, adviControl)
#' @noRd
.adviFitModel <- function(env) {
  stop("est=\"advi\" is not yet implemented", call. = FALSE)
}
