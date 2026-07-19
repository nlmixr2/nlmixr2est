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

# TRUE when the model carries a generalized (non-normal) / user-`ll()` endpoint.
# The npag Psi sums the inner per-observation llikObs, which for a non-normal
# endpoint is exactly the user's log-likelihood -- so the nonparametric objective
# is already correct; the residual-error (gamma) scaling is a no-op (r == 1).
#' @noRd
.npIsGeneralLik <- function(ui) {
  .dist <- tryCatch(ui$predDfFocei$distribution, error = function(e) NULL)
  if (is.null(.dist)) {
    .dist <- tryCatch(ui$predDf$distribution, error = function(e) NULL)
  }
  !is.null(.dist) && length(.dist) > 0L && any(as.character(.dist) != "norm")
}

# Auto initial-grid size when the user does not supply `points`.  Floors at 2028
# (the Pmetrics NPAG default, which covers a low-dimensional model well) and grows
# 512 points per support-point dimension (eta) beyond that, so a high-dimensional
# grid stays dense enough not to collapse.  Tested: theo (3 etas) -> 2028 matches
# Pmetrics; warfarin (8 etas) -> 4096 gives a non-degenerate PK/PD fit where a fixed
# 2028 collapsed.  The adaptive grid then refines from this seed.
#' @noRd
.npAutoPoints <- function(neta) {
  max(2028L, as.integer(512L * max(1L, as.integer(neta))))
}

# Validate the npag `dfScan` control: a single finite integer that is -1 (auto),
# 0 (skip the D(F) certificate), or a positive scan size.  Fails fast so a stray
# value (e.g. -5) cannot silently behave like auto.
#' @noRd
.npAssertDfScan <- function(dfScan) {
  if (length(dfScan) != 1L || is.na(dfScan) || !is.finite(dfScan) ||
        as.integer(dfScan) != dfScan || dfScan < -1L) {
    stop("'dfScan' must be a single integer: -1 (auto), 0 (skip), or a positive scan size",
         call. = FALSE)
  }
  as.integer(dfScan)
}

# Validate the npag/npb `cores` control: NULL (use the current rxode2 thread
# count, stored as NA) or a single positive integer thread count.
#' @noRd
.npAssertCores <- function(cores) {
  if (is.null(cores)) return(NA_integer_)
  if (length(cores) != 1L || is.na(cores) || !is.finite(cores) ||
        as.integer(cores) != cores || cores < 1L) {
    stop("'cores' must be NULL or a single positive integer number of threads",
         call. = FALSE)
  }
  as.integer(cores)
}

#' @noRd
.npEstCore <- function(env, est, muModel = NULL, ...) {
  .ui <- env$ui
  .what <- paste0(" for the estimation routine '", est, "'")
  # both npag and npb handle generalized (non-normal) / ll() endpoints: the
  # conditional likelihood sums the inner per-observation llikObs (correct for any
  # endpoint), which is all the npag grid and the npb Gibbs sweep need.
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
  # initial grid size: when the user does not supply `points`, scale it with the
  # number of support-point dimensions (etas).  A fixed grid (Pmetrics uses 2028)
  # covers a low-dimensional model well but grows sparse in high dimensions, where
  # the support can collapse; scale it so each added dimension keeps coverage.  npb
  # supplies its own K (truncation) so it never hits this auto path.
  .ctl$npPoints <-
    if (is.null(.ctl$points) || is.na(.ctl$points)) {
      .npAutoPoints(length(.box$names))
    } else as.integer(.ctl$points)
  .ctl$npCycles <- as.integer(if (is.null(.ctl$cycles)) 100L else .ctl$cycles)
  # global-optimality (D(F)) Sobol scan size (npag only): -1 auto
  # (max(2048, 2*npPoints)), 0 to skip the certificate, >0 for an explicit count.
  .ctl$npDfScan <- as.integer(if (is.null(.ctl$dfScan)) -1L else .ctl$dfScan)
  # gamma scales the residual variance r; a generalized (non-normal) endpoint has
  # r == 1, so the gamma warm-start is a no-op -- force it off there.
  .isGenLik <- .npIsGeneralLik(.ui)
  .ctl$npGammaOptimize <- !.isGenLik &&
    isTRUE(if (is.null(.ctl$gammaOptimize)) TRUE else .ctl$gammaOptimize)
  .ctl$npResidMode <- as.integer(
    if (is.null(.ctl$residOptimize)) 1L
    else switch(.ctl$residOptimize, none = 0L, alternate = 1L, final = 2L, 1L))
  # npb (stick-breaking Gibbs) knobs; npPoints doubles as the truncation level K
  .ctl$npAlpha <- as.numeric(if (is.null(.ctl$alpha)) 1.0 else .ctl$alpha)
  .ctl$npBurnin <- as.integer(if (is.null(.ctl$burnin)) 500L else .ctl$burnin)
  .ctl$npNsamp <- as.integer(if (is.null(.ctl$nsamp)) 500L else .ctl$nsamp)
  .ctl$npPropSd <- as.numeric(if (is.null(.ctl$propSd)) 0.2 else .ctl$propSd)
  .ctl$npSeed <- as.integer(if (is.null(.ctl$seed)) 42L else .ctl$seed)
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
  # honor npagControl(cores=)/npbControl(cores=): the parallel per-subject solves
  # in the kernels use the solve's thread count (op->cores), which rxode2 sizes
  # from getRxThreads() when the inner solve is built inside .npFamilyFit.  Set
  # the thread count for the fit (default: the current rxode2 threads) and restore
  # it afterwards; the fit results are independent of the thread count.
  .nCores <- if (is.null(.ctl$npCores) || is.na(.ctl$npCores)) {
    rxode2::getRxThreads()
  } else as.integer(.ctl$npCores)
  if (!is.na(.nCores) && .nCores >= 1L && rxode2::getRxThreads() != .nCores) {
    .oldThreads <- rxode2::getRxThreads()
    rxode2::setRxThreads(.nCores)
    on.exit(rxode2::setRxThreads(.oldThreads), add = TRUE)
  }
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
  # the nonparametric kernels do not compute the Monte-Carlo importance-sampling
  # ("imp") covariance; force it off (impmapControl's default) so the post-fit
  # recompute installs the FOCEI covariance instead
  .control$impCov <- FALSE
  .covMethodUser <- .control$covMethod  # restored on the fit env control below
  if (is.null(.covMethodUser) || identical(as.integer(.covMethodUser), 0L)) {
    # the inherited "imp" default maps to covMethod=0L for np; recompute the
    # analytic FOCEI covariance instead of reporting none
    .covMethodUser <- 2L
    .control$covType <- "analytic"
  }
  .control$maxOuterIterations <- 0L
  .control$covMethod <- 0L  # covariance is computed post-fit (.foceiRecomputeMuCov)
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
  # fixed thetas are held constant: exclude them from the mean-shift (they keep
  # their ini value instead of being moved to the support-point mean).
  .thOrd <- .th[order(.th$ntheta), , drop = FALSE]
  .thFixed <- which(!is.na(.thOrd$fix) & .thOrd$fix) - 1L   # 0-based fixed theta idx
  .keep <- !is.na(.muThetaIdx) & !is.na(.muEtaIdx) &
    !(.muThetaIdx %in% .covGroupTheta) & !(.muThetaIdx %in% .thFixed)
  .control$impMuThetaIdx <- .muThetaIdx[.keep]
  .control$impMuEtaIdx <- .muEtaIdx[.keep]
  .control$impThetaSensIdx <- integer(0)   # no sensitivity model for npag/npb
  .etaOrd <- .etaRows[order(.etaRows$neta1), ]
  .control$impOmegaFixedEta <- as.integer(which(isTRUE(.etaOrd$fix) | .etaOrd$fix) - 1L)
  # 0-based theta indices of the variance-scale residual parameters (add/prop/
  # lnorm/...).  The npag/npb assay-error multiplier (gamma) scales the residual
  # variance r; at finalization gamma is folded into these coefficients so the
  # reported parameter reflects the estimate.  Transform (boxCox/yeoJohnson) and
  # autocorrelation (ar) params are NOT variance scales and must not be folded.
  .errType <- as.character(.thOrd$err)
  .isFix <- !is.na(.thOrd$fix) & .thOrd$fix
  # variance-scale residual params (add/prop/lnorm/...) -- the ones the gamma
  # warm-start folds into.  Transform (boxCox/yeoJohnson), autocorrelation (ar)
  # and the power exponent (pw) are not variance scales.  Fixed ones are held.
  .errScale <- !is.na(.thOrd$err) &
    !(.errType %in% c("boxCox", "yeoJohnson", "ar", "pw")) & !.isFix
  .control$npResidScaleIdx <- as.integer(which(.errScale) - 1L)
  # ALL non-fixed residual (err-tagged) params are optimized by the residual step
  # (mixing distribution held fixed) against the extended-least-squares objective at
  # the posterior-mean etas.  kind selects the coordinate mapping: 1 = positive (SD),
  # 2 = correlation in (-1,1) (ar), 0 = free (lambda/exponent).
  .errOpt <- !is.na(.thOrd$err) & !.isFix
  .control$npResidOptIdx <- as.integer(which(.errOpt) - 1L)
  .optType <- .errType[.errOpt]
  .kind <- rep(1L, length(.optType))
  .kind[.optType %in% c("ar")] <- 2L
  .kind[.optType %in% c("boxCox", "yeoJohnson", "pw")] <- 0L
  .control$npResidOptKind <- as.integer(.kind)
  # per-residual-param endpoint (0-based, predDf row order) and proportional flag,
  # for the saem-style moment warm start.  The err row's `condition` names its
  # endpoint variable; a proportional term uses (err/f), everything else (additive,
  # lognormal, box-cox -- all additive on the transform-both-sides scale) uses err.
  .endVar <- tryCatch(as.character(ui$predDf$cond), error = function(e) character(0))
  .errCond <- as.character(.thOrd$condition[.errOpt])
  .errEnd <- match(.errCond, .endVar) - 1L
  .errEnd[is.na(.errEnd)] <- 0L
  .control$npResidOptEnd <- as.integer(.errEnd)
  .control$npResidOptProp <- as.integer(startsWith(.optType, "prop"))
  # per-endpoint compartment (predDf order) so the C++ side can map each observation's
  # cmt (rxode2 getIndCmt, the CMT time-varying covariate) to its endpoint index for the
  # per-endpoint moment warm start.  cmt values are distinct, not necessarily sequential.
  .control$npEndpointCmt <- tryCatch(as.integer(ui$predDf$cmt), error = function(e) integer(0))
  # ini-block bounds of the residual-opt params (for the bounded bobyqa step),
  # intersected with the parameter's natural range: an SD (kind 1) is >= 0 and the
  # continuous-AR correlation (kind 2) is in (0, 1).
  .lo <- as.numeric(.thOrd$lower[.errOpt])
  .hi <- as.numeric(.thOrd$upper[.errOpt])
  .lo[.kind == 1L] <- pmax(.lo[.kind == 1L], 0)
  .lo[.kind == 2L] <- pmax(.lo[.kind == 2L], 0)
  .hi[.kind == 2L] <- pmin(.hi[.kind == 2L], 0.999)
  .control$npResidOptLower <- .lo
  .control$npResidOptUpper <- .hi
  if (is.null(.control$npResidMode)) .control$npResidMode <- 1L
  .est <- if (exists("est", envir = env)) get("est", envir = env) else "npag"
  # mixture (mix()) proportions: npag estimates them via the in-cycle EM update; npb
  # samples them via a Dirichlet(alpha0 + component counts) Gibbs step.  Either way
  # the components are marginalized in the conditional likelihood.  Skip the update/
  # sampling (hold the ini proportions) only when every proportion is fixed.
  .mixNames <- ui$mixProbs
  .npMixOptimize <- FALSE
  if (length(.mixNames) > 0) {
    .mixRows <- .thOrd[.thOrd$name %in% .mixNames, , drop = FALSE]
    .npMixOptimize <- any(!(!is.na(.mixRows$fix) & .mixRows$fix))
  }
  .control$npMixOptimize <- isTRUE(.npMixOptimize)
  # freeze flag is now vestigial (the residual step re-solves at the posterior-mean
  # etas); kept for the C++ signature.  The err-tagged residual params feed only f/r.
  .control$npResidFreeze <- TRUE
  # A non-mu, non-fixed structural fixed-effect theta (err==NA, no eta) is neither a
  # grid dimension nor a residual/likelihood parameter.  With muExpand=TRUE it has
  # been injected as a pseudo-eta (so it is now mu-referenced and not seen here);
  # otherwise, optimize it as a "regressor" against the MARGINAL likelihood over the
  # whole support (a structural shift is not identified by the residual ELS step at
  # fixed etas), re-solving the ODE per candidate.  Kept SEPARATE from the residual
  # (ELS) set -- mixing the two objectives mis-identifies both.
  .isMuRef <- .thOrd$name %in% .mr$theta
  .isMixP <- .thOrd$name %in% .mixNames
  # A non-mu, non-fixed structural fixed-effect theta (err==NA, no eta) is optimized as
  # a "regressor".  This INCLUDES a mix() model's structural component parameters (e.g.
  # the per-component typical clearance): the residual/regressor objective is the exact
  # conditional likelihood, which for a mixture marginalizes over the components with the
  # held proportions (NONMEM7 eq 1.182), so moving a component parameter changes that
  # component's fit and the mixture likelihood -- the components are estimated, not held.
  # Only the mixture PROPORTIONS are excluded (they move in the proportion update step).
  .regress <- which(is.na(.thOrd$err) & !.isFix & !.isMuRef & !.isMixP)
  if (length(.regress) > 0L) {
    .control$npRegressIdx <- as.integer(.regress - 1L)
    .control$npRegressLower <- as.numeric(.thOrd$lower[.regress])
    .control$npRegressUpper <- as.numeric(.thOrd$upper[.regress])
    # residOptimize="none" (npResidMode 0) holds everything, regressors included --
    # do not override it here.
  } else {
    .control$npRegressIdx <- integer(0)
    .control$npRegressLower <- numeric(0)
    .control$npRegressUpper <- numeric(0)
  }
  # mu-expanded (injected) etas: finalization recovers each as a FIXED effect by
  # folding its support-mean into the paired theta and collapsing its random effect.
  # Pass the 0-based eta index (etaNames order) and theta index (thetaTrans/ntheta
  # order) of each injected pair to the C++ finalizer.
  .pairs <- nlmixr2global$npMuExpandPairs
  .injEta <- if (is.null(.pairs)) character(0) else .pairs$eta
  .injTh <- if (is.null(.pairs)) character(0) else .pairs$theta
  if (length(.injEta) > 0L) {
    .ei <- match(.injEta, .etaNames) - 1L
    .ti <- match(.injTh, .thNames) - 1L
    .ok <- !is.na(.ei) & !is.na(.ti)
    .control$npMuExpandEtaIdx <- as.integer(.ei[.ok])
    .control$npMuExpandThetaIdx <- as.integer(.ti[.ok])
  }
  assign("control", .control, envir = ui)
  .fit <- .foceiFamilyReturn(env, ui, ..., est = .est)
  .impRestoreCovMethod(.fit, .covMethodUser)
  .fit
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
  .np <- list(points = NA_integer_, cycles = 100L, gammaOptimize = TRUE,
              residOptimize = "alternate", muExpand = FALSE,
              gridWidth = 4, gridBounds = "auto", dfScan = -1L, npCores = NA_integer_,
              alpha = 1.0, burnin = 500L, nsamp = 500L, nchains = 1L,
              propSd = 0.2, seed = 42L)
  for (.n in names(.np)) if (!is.null(.in[[.n]])) .np[[.n]] <- .in[[.n]]
  if (is.list(.in)) {
    for (.n in c(names(.np), "est")) .in[[.n]] <- NULL   # strip npag/npb-only fields
  }
  .ctl <- getValidNlmixrCtl.impmap(list(.in))
  .ctl$est <- est
  .ctl$points <- as.integer(.np$points)
  .ctl$cycles <- as.integer(.np$cycles)
  .ctl$gammaOptimize <- isTRUE(.np$gammaOptimize)
  .ctl$residOptimize <- as.character(.np$residOptimize)
  .ctl$muExpand <- isTRUE(.np$muExpand)
  .ctl$gridWidth <- as.numeric(.np$gridWidth)
  .ctl$gridBounds <- as.character(.np$gridBounds)
  .ctl$dfScan <- as.integer(.np$dfScan)
  .ctl$npCores <- as.integer(.np$npCores)
  .ctl$alpha <- as.numeric(.np$alpha)
  .ctl$burnin <- as.integer(.np$burnin)
  .ctl$nsamp <- as.integer(.np$nsamp)
  .ctl$nchains <- as.integer(.np$nchains)
  .ctl$propSd <- as.numeric(.np$propSd)
  .ctl$seed <- as.integer(.np$seed)
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
