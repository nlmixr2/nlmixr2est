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
  .foceiInnerControl(
    control,
    likelihood = control$likelihood,
    literalFixRes = control$literalFixRes,
    eventSens = control$eventSens,
    indTolRelax = control$indTolRelax,
    stickyRecalcN = control$stickyRecalcN,
    # the analytic outer solve's own loosening: est="vae" reuses the
    # SAME inner call to choose the likelihood, so it gets the same
    # generalization and the same knobs rather than a parallel set
    outerMaxOdeRecalc = control$outerMaxOdeRecalc,
    outerOdeRecalcFactor = control$outerOdeRecalcFactor,
    outerStickyRecalcN = control$outerStickyRecalcN,
    # the per-subject FD step of the outer gradient's fallback: est="vae"
    # reaches the same code through foceiGradPooledDirect_, so it gets the
    # same knob rather than a parallel one
    fdIndividualStep = if (is.null(control$fdIndividualStep)) {
      TRUE
    } else {
      isTRUE(control$fdIndividualStep)
    },
    fdOutlierZ = if (is.null(control$fdOutlierZ)) {
      3.5
    } else {
      as.double(control$fdOutlierZ)
    },
    fdOutlierScale = if (is.null(control$fdOutlierScale)) {
      TRUE
    } else {
      isTRUE(control$fdOutlierScale)
    },
    fdRefine = if (is.null(control$fdRefine)) {
      "chartrand"
    } else {
      as.character(control$fdRefine)
    },
    fdChartrandAll = isTRUE(control$fdChartrandAll),
    fdOutlierAny = isTRUE(control$fdOutlierAny)
  )
}

#' Set up the FOCEi inner problem for `ui` at its current ini() estimates
#' (`.foceiInnerEnv()`), plus the VAE's own omega packing and, for
#' nonMuTheta="grad", the augmented outer model; then the C++ vaeInnerSetup_
#' (foceiSetup_ + updateTheta).  Keep the env alive until vaeInnerFree().
#' @noRd
.vaeInnerSetup <- function(ui, data, etaMat, control, est = "focei") {
  .env <- .foceiInnerEnv(rxode2::rxUiDecompress(ui), data, .vaeInnerFoceiControl(control), est, etaMat)
  .ui <- .env$ui
  ## A non-Gaussian endpoint has no eta-epsilon interaction term to carry: rx_pred_
  ## IS the log-density.  The focei flow pairs needOptimHess with interaction=0 for
  ## that reason (.foceiFitInternal); this entry must do the same, or the inner
  ## problem is set up for the FOCEi (f,R) kernel while the objective runs the
  ## exact-Hessian one.
  if (isTRUE(.env$control$needOptimHess)) {
    .env$control$interaction <- 0L
  }
  ## "sqrt"-xform rxInv on the model's DECLARED omega structure (diagonal plus
  ## any correlated blocks): the per-step C++ fast path (vaeInnerUpdatePar_)
  ## packs chol(Omega^-1) onto the omega block of the reduced par vector, using
  ## the 0-based position list stashed here (column-major upper-tri restricted
  ## to the structure -- rxSymInvCholCreate's parameter order)
  ## same repair ladder as focei: a 0 sitting inside a correlated block is not
  ## representable here, so it is filled and estimated instead of aborting the
  ## run with "theta has to have N elements" (#1079)
  ## the foceiOptEnv build already reported any repair
  .sic <- .foceiSymInvCholCreate(.ui$omega, "sqrt", NULL, warn = FALSE)
  .om <- .sic$mat
  .env$rxInv <- .sic$rxInv
  .selMat <- upper.tri(.om, diag = TRUE) & .om != 0
  diag(.selMat) <- TRUE
  .env$vaeOmegaSel <- which(.selMat, arr.ind = TRUE) - 1L
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
      ## Registering it on the model list is enough: the C++ pool registry sizes
      ## the pool for the largest peer and derives the inner override itself, so R
      ## no longer nominates a poolModel or computes innerNeq.  foceiSetup_ still
      ## aliases its THETA_1_/ETA_1_ spelling onto the THETA[1]/ETA[1] columns so
      ## rxSolve_ can bind it.
      .env$model$vaeOuter <- .am$augMod
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
  .cores <- tryCatch(
    {
      .c <- control$rxControl$cores
      if (is.null(.c) || is.na(.c) || .c < 1L) as.integer(rxode2::getRxThreads()) else as.integer(.c)
    },
    error = function(e) 1L
  )
  vaeInnerLik(as.matrix(etaMat), .cores, isTRUE(grad), isTRUE(preds))
}

#' Re-set up the inner problem at new population parameters (theta + omega)
#' without recompiling: reuses the env's compiled inner model and processed data,
#' rebuilds rxInv from the new omega, and re-runs foceiSetup_ + updateTheta.
#' @param env the setup env from .vaeInnerSetup
#' @param theta full theta vector (ntheta order): structural intercepts, error,
#'   covariate betas, mixture probs
#' @param omega random-effect variances: full matrix, or a vector taken as the
#'   diagonal (eta order)
#' @param etaMat starting etas [nsub, neta]
#' @noRd
.vaeInnerUpdate <- function(env, theta, omega, etaMat, diagXform = "sqrt") {
  env$thetaIni <- setNames(as.numeric(theta), paste0("THETA[", seq_along(theta), "]"))
  .om <- if (is.matrix(omega)) omega else diag(omega, length(omega))
  .nm <- env$etaNames
  if (!is.null(.nm) && length(.nm) == nrow(.om)) {
    dimnames(.om) <- list(.nm, .nm)
  }
  ## Reported once at setup, and this runs every VI step -- so no message, and
  ## no fallback either: only the block-zero fill (a 1e-10 correlation) may run
  ## here, a genuinely bad omega still errors rather than silently flooring.
  .sic <- .foceiSymInvCholCreate(.om, diagXform, NULL, warn = FALSE, fallback = FALSE)
  .om <- .sic$mat
  env$rxInv <- .sic$rxInv
  .selMat <- upper.tri(.om, diag = TRUE) & .om != 0
  diag(.selMat) <- TRUE
  env$vaeOmegaSel <- which(.selMat, arr.ind = TRUE) - 1L
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
.vaeElboStepInner <- function(
  params,
  prep,
  innerEnv,
  zPop,
  omega,
  a,
  alphaKL,
  eps,
  control,
  nMix = 1L,
  mixProb = 1,
  withGrad = TRUE
) {
  .cores <- tryCatch(
    {
      .c <- control$rxControl$cores
      if (is.null(.c) || is.na(.c) || .c < 1L) as.integer(rxode2::getRxThreads()) else as.integer(.c)
    },
    error = function(e) 1L
  )
  vaeElboStepCpp_(
    params,
    prep,
    zPop,
    if (is.matrix(omega)) omega else as.numeric(omega),
    as.numeric(a),
    as.numeric(alphaKL),
    as.matrix(eps),
    as.integer(nMix),
    as.numeric(mixProb),
    .cores,
    isTRUE(withGrad)
  )
}
