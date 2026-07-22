# vaeInner.R -- drive the FOCEi inner likelihood directly for the VAE. The inner
# problem is set up ONCE (foceiSetup_ via vaeInnerSetup_), then the per-subject
# (and per-mixture-component) objective/gradient are evaluated through
# likInner0/lpInner in a parallel OpenMP loop (vaeInnerLik) -- reusing the inner
# logic (mixtures, multiple endpoints, error structures, log-likelihood,
# censoring) without the nlmixr2 R interface.

#' A foceiControl carrying the VAE's chosen inner likelihood + solving options.
#' focei -> interaction=1; foce/focep -> interaction=0 (focep = FOCE+, R at the
#' live conditional eta); laplace -> the Laplace method.
#' @noRd
.vaeInnerFoceiControl <- function(control) {
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
  ## A non-Gaussian endpoint has no eta-epsilon interaction term to carry: rx_pred_
  ## IS the log-density.  The focei flow pairs needOptimHess with interaction=0 for
  ## that reason (.foceiFitInternal); this entry must do the same, or the inner
  ## problem is set up for the FOCEi (f,R) kernel while the objective runs the
  ## exact-Hessian one.
  if (isTRUE(.env$control$needOptimHess)) .env$control$interaction <- 0L
  ## AGQ off
  .env$aqn <- 0L; .env$qx <- double(0); .env$qw <- double(0); .env$qfirst <- FALSE
  .env$nAGQ <- 0L; .env$aqLow <- -Inf; .env$aqHi <- Inf; .env$nEstOmega <- 0L
  .env$etaMat <- etaMat
  ## force a diagonal "sqrt"-xform rxInv (the training parameterization): the
  ## per-step C++ fast path (vaeInnerUpdatePar_) maps eta variances onto the
  ## omega block of the reduced par vector, which requires omegan == neta
  .om <- .ui$omega
  .om <- diag(diag(.om), nrow(.om))
  dimnames(.om) <- dimnames(.ui$omega)
  .env$rxInv <- rxode2::rxSymInvCholCreate(mat = .om, diag.xform = "sqrt")
  ## nonMuTheta="grad": the augmented outer-gradient model is solved in the SHARED
  ## pool, so it must SIZE that pool -- it is the larger structure (26 states / 29
  ## lhs vs 6 / 6 on a one-compartment fit).  The inner MAP then runs under
  ## ind->neqOverride, exactly as est="impmap" does with its theta-sens model.
  ## Nothing is freed by the M-step, so no solve-arg stash is needed.
  if (identical(control$nonMuTheta, "grad")) {
    ## .ui$control was replaced with the DERIVED focei control above, so
    ## .analyticGradCaller (which rxUiGet.foceiOuter consults) would resolve to NA.
    ## Re-mark it before asking for the augmented model.
    .fcg <- .ui$control
    .fcg$nonMuTheta <- "grad"
    assign("control", .fcg, envir = .ui)
    .am <- tryCatch(.ui$foceiOuter, error = function(e) NULL)
    if (!is.null(.am) && inherits(.am$augMod, "rxode2") && !is.null(.env$model)) {
      .env$model$vaeOuter <- .am$augMod
      ## The augmented model SIZES the shared pool (26 states / 29 lhs vs the
      ## inner model's 6 / 6); the inner MAP then runs under ind->neqOverride.
      ## foceiSetup_ aliases its THETA_1_/ETA_1_ spelling onto the THETA[1]/ETA[1]
      ## columns so rxSolve_ can bind it.
      .env$poolModel <- .am$augMod
      .env$innerNeq <- length(rxode2::rxModelVars(.env$model$inner)$state)
    }
  }
  vaeInnerSetup_(.env)
  .env
}

#' Free the inner-problem state set up by .vaeInnerSetup.
#' @noRd
.vaeInnerFree <- function() invisible(vaeInnerFree_())

#' Evaluate the inner objective (and optionally the eta-gradient) at `etaMat`
#' (rows = ids: nSub, or nSub*nMix for mixtures) through the parallel C++ driver.
#' @noRd
.vaeInnerEval <- function(etaMat, control, grad = FALSE, preds = FALSE) {
  .cores <- tryCatch({
    .c <- control$rxControl$cores
    if (is.null(.c) || is.na(.c) || .c < 1L) as.integer(rxode2::getRxThreads()) else as.integer(.c)
  }, error = function(e) 1L)
  vaeInnerLik(as.matrix(etaMat), .cores, isTRUE(grad), isTRUE(preds))
}

#' Re-set up the inner problem at new population parameters (theta + omega)
#' without recompiling: reuses the env's compiled inner model and processed data,
#' rebuilds rxInv from the new omega, and re-runs foceiSetup_ + updateTheta.
#' @param env the setup env from .vaeInnerSetup
#' @param theta full theta vector (ntheta order): structural intercepts, error,
#'   covariate betas, mixture probs
#' @param omega diagonal random-effect variances (eta order)
#' @param etaMat starting etas [nsub, neta]
#' @noRd
.vaeInnerUpdate <- function(env, theta, omega, etaMat, diagXform = "sqrt") {
  env$thetaIni <- setNames(as.numeric(theta), paste0("THETA[", seq_along(theta), "]"))
  .om <- diag(omega, length(omega))
  .nm <- env$etaNames
  if (!is.null(.nm) && length(.nm) == nrow(.om)) dimnames(.om) <- list(.nm, .nm)
  env$rxInv <- rxode2::rxSymInvCholCreate(mat = .om, diag.xform = diagXform)
  env$etaMat <- etaMat
  vaeInnerSetup_(env)
  invisible(env)
}

#' One ELBO evaluation using the FOCEi inner likelihood -- a thin R interface to
#' the C++ core (`vaeElboStepCpp_`) that the training loop (`vaeTrainCpp_`) also
#' calls directly. likInner0(eta) = p(x|z) + p(z); its eta-gradient is the encoder
#' upstream gZ = lp - Omega^-1 eta + alphaKL*(z - z_pop)/Omega, gLogSigma =
#' -alphaKL. Mixtures (nMix>1) evaluate nSub*nMix ids and combine with
#' -2 logsumexp over mixProb. The inner problem must already be set up
#' (`.vaeInnerSetup`); `innerEnv` is accepted for signature compatibility but the
#' C++ core reads the active op_focei allocation that setup created.
#' @noRd
.vaeElboStepInner <- function(params, prep, innerEnv, zPop, omega, a, alphaKL, eps,
                              control, nMix = 1L, mixProb = 1, withGrad = TRUE) {
  .cores <- tryCatch({
    .c <- control$rxControl$cores
    if (is.null(.c) || is.na(.c) || .c < 1L) as.integer(rxode2::getRxThreads()) else as.integer(.c)
  }, error = function(e) 1L)
  vaeElboStepCpp_(params, prep, zPop, as.numeric(omega), as.numeric(a),
                  as.numeric(alphaKL), as.matrix(eps), as.integer(nMix),
                  as.numeric(mixProb), .cores, isTRUE(withGrad))
}
