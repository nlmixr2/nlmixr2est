# Methods for newuoa residual optimization
.saemResidF <- function(x) {
  # Exports the C residual functions in saem
  .Call(`_saemResidF`, x)
}

.saemOpt1 <- function(p1) {
  .opt <- stats::nlm(.saemResidF, p1)
  .opt$estimate
}

.newuoa <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .ctl <- control
  if (is.null(.ctl$npt)) .ctl$npt <- length(par) * 2 + 1
  if (is.null(.ctl$rhobeg)) .ctl$rhobeg <- 0.2
  if (is.null(.ctl$rhoend)) .ctl$rhoend <- 1e-4
  .ctl$iprint <- 0L
  .ctl <- .ctl[names(.ctl) %in% c("npt", "rhobeg", "rhoend", "iprint", "maxfun")]
  .ret <- try(minqa::newuoa(par, fn,
    control = .ctl
    ))
  if (inherits(.ret, "try-error")) {
    .ret <- list()
    .ret$x <- rep(NA_real_, length(.ret$par))
    .ret$message <- "try-error"
    .ret$convergence <- -42L
    .ret$value <- NA_real_
    return(.ret)
  }
  .ret$x <- .ret$par
  .ret$message <- .ret$msg
  .ret$convergence <- .ret$ierr
  .ret$value <- .ret$fval
  return(.ret)
}

.saemCheckCfg <- function(cfg) {
  checkmate::assertIntegerish(cfg$itmax, lower=1, len=1, .var.name="saem.cfg$itmax")
  checkmate::assertNumeric(cfg$tol, lower=0, len=1, .var.name="saem.cfg$tol")
  # currently supports 1=nealder meade and 2=newoua
  checkmate::assertIntegerish(cfg$type, lower=1, upper=2, len=1, .var.name="saem.cfg$type")
  checkmate::assertNumeric(cfg$lambdaRange, lower=0, len=1, .var.name="saem.cfg$lambdaRange")
  checkmate::assertIntegerish(cfg$maxOdeRecalc, lower=0, len=1, .var.name="saem.cfg$maxOdeRecalc")
  checkmate::assertNumeric(cfg$odeRecalcFactor, lower=0, len=1, .var.name="saem.cfg$odeRecalcFactor")
  # nmc number of mc interations
  checkmate::assertIntegerish(cfg$nmc, lower=0, len=1, .var.name="saem.cfg$nmc")
  .nmc <- cfg$nmc
  # nu is the number of selection for each probability type.
  checkmate::assertIntegerish(cfg$nu, lower=0, len=3, .var.name="saem.cfg$nu")
  # Overall number of iterations
  checkmate::assertIntegerish(cfg$niter, lower=0,  len=1, .var.name="saem.cfg$niter")
  # Number of iterations where the correlation is ignored
  checkmate::assertIntegerish(cfg$nb_correl, lower=0,  len=1, .var.name="saem.cfg$nb_correl")
  # Number of iteractionons for phi0
  checkmate::assertIntegerish(cfg$niter_phi0, lower=0, len=1, .var.name="saem.cfg$niter_phi0")
  checkmate::assertNumeric(cfg$coef_phi0, lower=0, len=1, .var.name="saem.cfg$coef_phi0")
  # Coefficient for SA step and number of SA iterations
  checkmate::assertIntegerish(cfg$nb_sa, lower=0, len=1, .var.name="saem.cfg$nb_sa")
  checkmate::assertNumeric(cfg$coef_sa, lower=0, len=1, .var.name="saem.cfg$coef_sa")
  #
  checkmate::assertNumeric(cfg$rmcmc, lower=0, len=1, .var.name="saem.cfg$rmcmc")
  # Length of pas, pash = niter, both factors for reduction in saem
  checkmate::assertNumeric(cfg$pas, lower=0, .var.name="saem.cfg$pas")
  checkmate::assertNumeric(cfg$pash, lower=0, .var.name="saem.cfg$pash")
  # Number of parameters with thetas and etas (nphi1) and simply thetas (nphi0)
  checkmate::assertIntegerish(cfg$nphi1 , lower=0, .var.name="saem.cfg$nphi1")
  checkmate::assertIntegerish(cfg$nphi0 , lower=0, .var.name="saem.cfg$nphi0")
  .nphi1 <- cfg$nphi1
  .nphi0 <- cfg$nphi0
  .nphi <- .nphi0 + .nphi1
  # minimum value of a parameter
  checkmate::assertNumeric(cfg$minv, lower=0, len=.nphi, .var.name="saem.cfg$minv")
  # N is the number of IDs
  checkmate::assertIntegerish(cfg$N, lower=0, len=1, .var.name="saem.cfg$N")
  .N <- cfg$N
  # Total number of items in the dataset
  checkmate::assertIntegerish(cfg$ntotal, lower=0, len=1, .var.name="saem.cfg$ntotal")
  .ntotal <- cfg$ntotal
  # observed
  checkmate::assertNumeric(cfg$y, len=.ntotal, .var.name="saem.cfg$y")

  # event table matrix
  checkmate::assertMatrix(cfg$evt, mode="numeric", .var.name="saem.cfg$evt")
  # phi matrix
  checkmate::assertMatrix(cfg$phiM, mode="numeric", ncols=.nphi, .var.name="saem.cfg$phiM")

  checkmate::assertIntegerish(cfg$mlen,  lower=1, len=1, .var.name="saem.cfg$mlen")

  # maximum number of measurments for an indiviaul
  .mlen <- cfg$mlen

  checkmate::assertIntegerish(cfg$indio, min.len=1, .var.name="saem.cfg$indio")

  # covstruct and Mcovariables
  checkmate::assertMatrix(cfg$covstruct1, mode="numeric", .var.name="saem.cfg$covstruct1")
  checkmate::assertMatrix(cfg$Mcovariables, mode="numeric", .var.name="saem.cfg$Mcovariables")

  # indicator colum for phi1 and phi0
  checkmate::assertIntegerish(cfg$i1, .var.name="saem.cfg$i1")

  checkmate::assertMatrix(cfg$Gamma2_phi1, mode="numeric", .var.name="saem.cfg$Gamma2_phi1")
  checkmate::assertMatrix(cfg$Gamma2_phi1fixedIx, mode="integerish", .var.name="saem.cfg$Gamma2_phi1fixedIx")

  .Gamma2_phi1fixed <- cfg$Gamma2_phi1fixed
  # phi1fixed indicator
  checkmate::assertIntegerish(cfg$Gamma2_phi1fixed,lower=0, upper=1, len=1, .var.name="saem.cfg$Gamma2_phi1fixed")
  if (.Gamma2_phi1fixed == 1) {
    checkmate::assertMatrix(cfg$Gamma2_phi1fixedValues, mode="numeric", .var.name="saem.cfg$Gamma2_phi1fixedValues")
  }
  checkmate::assertMatrix(cfg$mprior_phi1, mode="numeric", .var.name="saem.cfg$mprior_phi1")

  checkmate::assertMatrix(cfg$COV1, mode="numeric", .var.name="saem.cfg$COV1")
  checkmate::assertMatrix(cfg$LCOV1, mode="numeric", .var.name="saem.cfg$LCOV1")
  checkmate::assertMatrix(cfg$COV21, mode="numeric", .var.name="saem.cfg$COV21")
  checkmate::assertMatrix(cfg$MCOV1, mode="numeric", .var.name="saem.cfg$MCOV1")

  checkmate::assertIntegerish(cfg$jcov1, .var.name="saem.cfg$jcov1")
  checkmate::assertIntegerish(cfg$ind_cov1, .var.name="saem.cfg$ind_cov1")

  checkmate::assertMatrix(cfg$statphi11, mode="numeric", .var.name="cfg$statphi11")
  checkmate::assertMatrix(cfg$statphi12, mode="numeric", .var.name="cfg$statphi12")
  if (.nphi0 > 0) {
    checkmate::assertIntegerish(cfg$i0, .var.name="saem.cfg$i0")
    checkmate::assertMatrix(cfg$Gamma2_phi0, mode="numeric", .var.name="saem.cfg$Gamma2_phi0")
    checkmate::assertMatrix(cfg$mprior_phi0, mode="numeric", .var.name="saem.cfg$mprior_phi0")
    checkmate::assertMatrix(cfg$COV0, mode="numeric", .var.name="saem.cfg$COV0")
    checkmate::assertMatrix(cfg$LCOV0, mode="numeric", .var.name="saem.cfg$LCOV0")
    checkmate::assertMatrix(cfg$COV20, mode="numeric", .var.name="saem.cfg$COV20")
    checkmate::assertMatrix(cfg$MCOV0, mode="numeric", .var.name="saem.cfg$MCOV0")
    checkmate::assertIntegerish(cfg$jcov0, .var.name="saem.cfg$jcov0")
    checkmate::assertIntegerish(cfg$ind_cov0, .var.name="saem.cfg$ind_cov0")
    checkmate::assertMatrix(cfg$statphi01, mode="numeric", .var.name="cfg$statphi01")
    checkmate::assertMatrix(cfg$statphi02, mode="numeric", .var.name="cfg$statphi02")
  }

  checkmate::assertIntegerish(cfg$nlambda1, len=1, .var.name="saem.cfg$nlambda1")
  checkmate::assertIntegerish(cfg$nlambda0, len=1, .var.name="saem.cfg$nlambda0")
  checkmate::assertIntegerish(cfg$ilambda1, .var.name="saem.cfg$ilambda1")
  checkmate::assertIntegerish(cfg$ilambda0, .var.name="saem.cfg$ilambda0")

  checkmate::assertIntegerish(cfg$nendpnt, len=1, .var.name="saem.cfg$nendpnt")
  .nendpnt <- cfg$nendpnt
  checkmate::assertIntegerish(cfg$ix_sorting, .var.name="ix_sorting")
  checkmate::assertNumeric(cfg$ys, .var.name="saem.cfg$ys")
  checkmate::assertIntegerish(cfg$y_offset, .var.name="saem.cfg$y_offset")
  # The should match the number of endpoints
  checkmate::assertIntegerish(cfg$res.mod, len=.nendpnt, .var.name="saem.cfg$res.mod")
  checkmate::assertNumeric(cfg$ares, len=.nendpnt, .var.name="saem.cfg$ares")
  checkmate::assertNumeric(cfg$bres, len=.nendpnt, .var.name="saem.cfg$bres")
  checkmate::assertNumeric(cfg$cres, len=.nendpnt, .var.name="saem.cfg$cres")
  checkmate::assertNumeric(cfg$lres, len=.nendpnt, .var.name="saem.cfg$lres")
  checkmate::assertIntegerish(cfg$yj, lower=0, len=.nendpnt, .var.name="saem.cfg$yj")
  checkmate::assertIntegerish(cfg$propT, lower=0, upper=1, len=.nendpnt, .var.name="saem.cfg$yj")
  checkmate::assertNumeric(cfg$lambda, len=.nendpnt, .var.name="cfg$lambda")
  checkmate::assertNumeric(cfg$low, len=.nendpnt, .var.name="cfg$low")
  checkmate::assertNumeric(cfg$hi, len=.nendpnt, .var.name="cfg$hi")
  checkmate::assertNumeric(cfg$addProp, len=.nendpnt, .var.name="cfg$addProp")
  checkmate::assertIntegerish(cfg$ix_endpnt, .var.name="cfg$ix_endpnt")
  checkmate::assertMatrix(cfg$ix_idM, mode="integerish", .var.name="cfg$ix_idM")
  checkmate::assertIntegerish(cfg$res_offset, .var.name="cfg$res_offset")

  checkmate::assertIntegerish(cfg$nb_fixOmega, len=1, lower=0, .var.name="cfg$nb_fixOmega")
  checkmate::assertIntegerish(cfg$nb_fixResid, len=1, lower=0, .var.name="cfg$nb_fixResid")

  checkmate::assertNumeric(cfg$resValue, .var.name="cfg$resValue")
  checkmate::assertIntegerish(cfg$resFixed, lower=0, upper=1, .var.name="cfg$resFixed")
}

#' Fit a UI model with saem
#'
#' @param ui rxode2 ui
#' @param data nlmixr data
#' @param timeVaryingCovariates Time varying covarites in the data
#' @return lower level saem fit
#' @author Matthew L. Fidler
#' @noRd
.saemFitModel <- function(ui, data, timeVaryingCovariates=character(0)) {
  # Stage the time-varying/non-time-varying mu-ref covariate split so
  # $saemModel0 collapses to the phi + timeVaryingCovariate*beta model (shared
  # with the mu-referenced focei family, see .nlmixrSetMuRefTimeVarying).
  .nlmixrSetMuRefTimeVarying(ui, timeVaryingCovariates)
  on.exit(.nlmixrRmMuRefTimeVarying(ui))
  # Building the saem model list does the symengine translation and rxode2
  # compilation -- timed as "configure" (mapped to "setup") so it is not
  # silently absorbed into the "other" bucket.
  .model <- nlmixrWithTiming("configure", {
    ui$saemModelList
  })
  .inits <- ui$saemInit
  .rxControl <- rxode2::rxGetControl(ui, "rxControl", rxode2::rxControl())
  ## Delay differential equation models need a dense-output solver so delay()
  ## history is RECORDED (dense=TRUE) and interpolated; the SAEM default
  ## (liblsoda/lsoda) does neither, which mis-evaluates delay() throughout the fit
  ## and yields a non-finite covariance linearization.  Mirror rxode2::rxSolve()'s
  ## hasDelay enforcement here so the SAEM solve and the covariance dopred both use
  ## the dense dop853 path.
  if (isTRUE(rxode2::rxModelVars(attr(.model$saem_mod, "rx"))$flags[["hasDelay"]] == 1L)) {
    .rxControl$method <- 0L  # dop853 (dense; no analytic Jacobian required)
    .rxControl$stiff2 <- 0L
    .rxControl$dense <- TRUE # record dense history for delay() interpolation
  }
  .ue <- .uninformativeEtas(ui,
                            handleUninformativeEtas=rxode2::rxGetControl(ui, "handleUninformativeEtas", TRUE),
                            data=data,
                            attr(.model$saem_mod, "rx"),
                            rxControl=.rxControl)
  .seed <- rxode2::rxGetControl(ui, "seed", 99)
  # Run the whole fit -- the config-time R draws (rnorm) AND the C++ SAEM loop /
  # f-SAEM IMH kernel, which both draw through rxode2's threefry engine -- inside
  # rxWithSeed.  It sets BOTH R's RNG and the rxode2 engine seed and restores them
  # afterward, so the first fit of a session is properly seeded and fits never
  # contaminate each other's seed state (replaces set.seed + low-level seedEng()).
  rxode2::rxWithSeed(.seed, rxseed = .seed, {
  .cfg <- nlmixrWithTiming("configure", {
    .cfg <- .configsaem(model=.model,
                        data=data,
                        inits=.inits,
                        mcmc=rxode2::rxGetControl(ui, "mcmc",
                                                  list(niter = c(200, 300),
                                                       nmc = 3, nu = c(2, 2, 2))),
                        rxControl=.rxControl,
                        distribution=if (any(ui$predDf$distribution == "LL")) "general" else "normal",
                        fixedOmega=ui$saemModelOmegaFixed,
                        fixedOmegaValues=ui$saemModelOmegaFixedValues,
                        parHistThetaKeep=ui$saemParHistThetaKeep,
                        parHistOmegaKeep=ui$saemParHistOmegaKeep,
                        parHistOmegaOffPairs={
                          .oi <- ui$saemParHistOmegaOffInfo
                          if (is.null(.oi)) matrix(integer(0), ncol=2L) else .oi$pairs
                        },
                        seed=rxode2::rxGetControl(ui, "seed", 99),
                        DEBUG=rxode2::rxGetControl(ui, "DEBUG", 0),
                        tol=rxode2::rxGetControl(ui, "tol", 1e-6),
                        itmax=rxode2::rxGetControl(ui, "itmax", 30),
                        type=rxode2::rxGetControl(ui, "type", "newuoa"),
                        lambdaRange=rxode2::rxGetControl(ui, "lambdaRange", 3),
                        powRange=rxode2::rxGetControl(ui, "powRange", 10),
                        odeRecalcFactor=rxode2::rxGetControl(ui, "odeRecalcFactor", 10^0.5),
                        maxOdeRecalc=rxode2::rxGetControl(ui, "maxOdeRecalc", 10^0.5),
                        indTolRelax=rxode2::rxGetControl(ui, "indTolRelax", TRUE),
                        nSaCov=if (identical(rxode2::rxGetControl(ui, "covMethod", "linFim"), "sa"))
                                 as.integer(rxode2::rxGetControl(ui, "nSaCov", 500L)) else 0L,
                        nres=ui$saemModNumEst,
                        perSa=rxode2::rxGetControl(ui, "perSa", 0.75),
                        perNoCor=rxode2::rxGetControl(ui, "perNoCor", 0.75),
                        perFixOmega=rxode2::rxGetControl(ui, "perFixOmega", 0.1),
                        perFixResid=rxode2::rxGetControl(ui, "perFixResid", 0.1),
                        resFixed=as.integer(ui$saemResFixed),
                        ue=.ue,
                        mixProb=ui$saemMixProb,
                        mixProbMethod=rxode2::rxGetControl(ui, "mixProbMethod", "regress"),
                        mixProbStepExp=rxode2::rxGetControl(ui, "mixProbStepExp", 1),
                        mixProbPriorN=rxode2::rxGetControl(ui, "mixProbPriorN", 20),
                        mixSampleMethod=rxode2::rxGetControl(ui, "mixSampleMethod", "parallel"),
                        omegaShare=ui$saemOmegaShare,
                        omegaShareSubpop=ui$saemOmegaShareSubpop,
                        fast=rxode2::rxGetControl(ui, "fast", FALSE),
                        fastIter=rxode2::rxGetControl(ui, "fastIter", 20L),
                        # a general log-likelihood endpoint has no normal do_mcmc
                        # fallback, so the fast kernel must run every iteration
                        fastKernel=if (.fsaemGeneralLik(ui)) "throughout"
                                   else rxode2::rxGetControl(ui, "fastKernel", "firstN"),
                        fastCov=rxode2::rxGetControl(ui, "fastCov", "auto"),
                        fastLik=rxode2::rxGetControl(ui, "fastLik", "focei"))
    .cfg$nonMuTheta <- rxode2::rxGetControl(ui, "nonMuTheta", "regress")
    # integer gate the SAEM C++ reads: when 1, non-mu (phi0) thetas are
    # estimated by the bounded direct optimizer (bounds from phi0Lower/Upper)
    # for normal models too, not just general-likelihood.
    .cfg$nonMuThetaRegress <- as.integer(identical(.cfg$nonMuTheta, "regress"))
    # warm-start residual params from observed per-endpoint moments (npag-style)
    .cfg$residWarmStart <- as.integer(rxode2::rxGetControl(ui, "residWarmStart", TRUE))
    # mixProbMethod="regress": fix per-subject mixture membership (hard classify
    # once) instead of the soft-EM responsibility step.
    .cfg$mixProbRegress <- as.integer(identical(
      rxode2::rxGetControl(ui, "mixProbMethod", "regress"), "regress"))
    .cfg$cres <- ui$saemCres
    .cfg$yj <- ui$saemYj
    .cfg$lres <- ui$saemLres
    .cfg$low <- ui$saemLow
    .cfg$hi <- ui$saemHi
    .cfg$propT <- ui$saemPropT
    .cfg$addProp <- ui$saemAddProp
    # a general log-likelihood endpoint has no residual params, so these are
    # empty; coerce NULL to the typed empty the config checks expect
    .cfg$resValue <- as.numeric(ui$saemResValue)
    # Iteration-print flows through the shared scaleApplyIterPrintControl
    # (src/scale.h); U auto-skips since Plambda has no optimizer scaling,
    # leaving # and X. xform (xPar/probitIdx/bounds) drives the X row's
    # back-transform via scaleAttachXform, same as every other estimator.
    .cfg$parHistNames <- as.character(ui$saemParHistNames)
    .cfg$xform        <- .iterPrintXParFromUi(ui, .cfg$parHistNames)
    .cfg$iterPrintControl <- rxode2::rxGetControl(ui, "iterPrintControl",
                                                  iterPrintControl())
    # The general-likelihood phi0 step is optimized with the bounded bobyqa
    # (.boundedResidOpt) and the ODE states frozen.  Provide the phi0 (fixed-
    # effect-only) parameter bounds, in phi0 (i0) column order, from the theta
    # iniDf so the optimization stays in a valid region.
    if (!is.null(.cfg$nphi0) && .cfg$nphi0 > 0L) {
      .pars <- ui$saemParamsToEstimate
      .phi0Names <- .pars[.cfg$i0]
      .lo <- ui$iniDf$lower[match(.phi0Names, ui$iniDf$name)]
      .hi <- ui$iniDf$upper[match(.phi0Names, ui$iniDf$name)]
      .cfg$phi0Lower <- ifelse(is.na(.lo), -Inf, .lo)
      .cfg$phi0Upper <- ifelse(is.na(.hi), Inf, .hi)
    }
    if (isTRUE(rxode2::rxGetControl(ui, "fast", FALSE))) {
      .cfg <- .fsaemInstallStep(ui, data, .rxControl, .cfg)
    }
    .saemCheckCfg(.cfg)
    .cfg
  })
  .saemRes <- nlmixrWithTiming("saem", {
    .model$saem_mod(.cfg)
  })
  # f-SAEM sets up the FOCEi inner (op_focei globals + a shared solve); tear it
  # down so it does not leak into a later fit's solve state (reproducibility).
  if (isTRUE(rxode2::rxGetControl(ui, "fast", FALSE))) {
    try(vaeInnerFree_(), silent = TRUE)
  }
  .saemRes
  })
}
#' Get the saem control statement and install it into the ui
#'
#' @param env Environment with ui in it
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.saemFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- saemControl()
  }
  if (!inherits(.control, "saemControl")){
    .control <- do.call(nlmixr2est::saemControl, .control)
  }
  assign("control", .control, envir=.ui)
}
#' Get SAEM theta
#'
#' @param env Environment that has ui and saem in it
#' @return Nothining, environment is assigned a theta
#' @author Matthew L. Fidler
#' @noRd
.getSaemTheta <- function(env) {
  .ui <- env$ui
  .saem <- env$saem
  .iniDf <- .ui$iniDf
  .predDf <- .ui$predDf
  .thetaNames <- .iniDf[!is.na(.iniDf$ntheta), "name"]
  .theta <- setNames(rep(NA_real_, length(.thetaNames)), .thetaNames)
  .saemThetaNames <- .ui$saemParamsToEstimate
  .thetaSaem <- setNames(as.vector(fixed.effects(env$saem)), .saemThetaNames)
  .resMat <- .saem$resMat
  .varSpec <- FALSE
  for (n in .thetaNames) {
    if (n %in% .saemThetaNames) {
      .theta[n] <- .thetaSaem[n]
    }
  }
  .hasVariance <- any(names(.predDf) == "variance")
  for (i in seq_along(.predDf$cond)) {
    .x <- paste(.predDf$cond[i])
    .tmp <- .iniDf[which(.iniDf$condition == .x), ]
    .w <- which(vapply(.tmp$err,
                       function(x) any(x == c("prop", "propT", "pow", "powT")),
                       logical(1),
                       USE.NAMES=FALSE))
    if (length(.w) == 1) {
      .var <- FALSE
      if (.hasVariance) {
        .var <- .predDf$variance[i]
      }
      if (.var) {
        .theta[paste(.tmp$name[.w])] <- .resMat[i, 2]^2
        .varSpec <- TRUE
      } else {
        .theta[paste(.tmp$name[.w])] <- .resMat[i, 2]
      }
    }
    .w <- which(vapply(.tmp$err,
                       function(x) any(x == c("pow2", "powT2")),
                       logical(1),
                       USE.NAMES=FALSE))
    if (length(.w) == 1) {
      .theta[paste(.tmp$name[.w])] <- .resMat[i, 3]
    }
    .w <- which(vapply(seq_along(.tmp$err),
                       function(x) {
                         .x <- .tmp$err[x]
                         if (any(.x == c(
                           "add", "norm", "dnorm", "lnorm", "dlnorm",
                           "dlogn", "logn"))) {
                           if (!is.na(.tmp$est[x])) {
                             return(TRUE)
                           }
                         }
                         FALSE
                       },
                       logical(1),
                       USE.NAMES=FALSE))
    if (length(.w) == 1) {
      .var <- FALSE
      if (.hasVariance) {
        .var <- .predDf$variance[i]
      }
      if (.var) {
        .theta[paste(.tmp$name[.w])] <- .resMat[i, 1]^2
        .varSpec <- TRUE
      } else {
        .theta[paste(.tmp$name[.w])] <- .resMat[i, 1]
      }
    }
    .w <- which(vapply(.tmp$err, function(x) {
      any(x == c("boxCox", "yeoJohnson"))
    }, logical(1), USE.NAMES=FALSE))
    if (length(.w) == 1) {
      .theta[paste(.tmp$name[.w])] <- .resMat[i, 4]
    }
    # AR(1) autocorrelation estimated by the whitened M-step in src/saem.cpp
    .w <- which(vapply(.tmp$err, function(x) any(x == "ar"),
                       logical(1), USE.NAMES=FALSE))
    if (length(.w) == 1 && !is.null(.saem$arCor)) {
      .theta[paste(.tmp$name[.w])] <- .saem$arCor[i]
    }
  }
  if (length(.ui$mixProbs) > 0 && !is.null(.saem$mixProb)) {
    .estMix <- .saem$mixProb[seq_along(.ui$mixProbs)]
    .estMixClamped <- pmax(1e-6, pmin(1 - 1e-6, .estMix))
    .sumP <- sum(.estMixClamped)
    if (.sumP >= 1.0) {
      .estMixClamped <- .estMixClamped / (.sumP + 1e-6)
    }
    # Check collapse against the raw (pre-clamp) estimate, not the
    # clamp/raw difference: a component collapsed near 0 (e.g. 1e-8)
    # changes by only ~1e-6 in absolute terms when clamped to 1e-6, so
    # comparing the clamp delta against a tolerance misses the most
    # common collapse case. Rescaling (sum >= 1) is checked separately.
    .collapsed <- any(.estMix < 1e-3 | .estMix > 1 - 1e-3)
    .rescaled <- .sumP >= 1.0
    if (.collapsed || .rescaled) {
      warning("one or more estimated mixture probabilities collapsed toward 0/1 or required ",
              "rescaling; this may indicate the mixture components are not well identified",
              call. = FALSE)
    }
    .theta[.ui$mixProbs] <- .estMixClamped
  }

  env$fullTheta <- .theta
  if (.varSpec) {
    .minfo("variance residual estimates transformed from standard deviation")
    warning("variance residual estimates transformed from standard deviation",
            call.=FALSE)
  }
  invisible()
}
#' Get SAEM omega
#'
#' @param env Environment that has ui and saem in it
#' @return Nothing, environment is assigned the omega
#' @author Matthew L. Fidler
#' @noRd
.getSaemOmega <- function(env) {
  ## Reorder based on translation
  .saem <- env$saem
  .ui <- env$ui
  .etaTrans <- .ui$saemOmegaTrans
  ## saem eta ->  ui eta
  .df <- .ui$iniDf
  .eta <- .df[!is.na(.df$neta1), ]
  .etaNames <- .eta[.eta$neta1 == .eta$neta2, "name"]
  .neta <- length(.etaNames)
  .len <- length(.etaNames)
  .ome <- matrix(rep(0, .len * .len), .len, .len, dimnames=list(.etaNames, .etaNames))
  # Gamma2_phi1Report is the reporting-only pooled BSV for split ETAs; falls
  # back to Gamma2_phi1 for older cached fits without the field.
  .curOme <- if (!is.null(.saem$Gamma2_phi1Report)) .saem$Gamma2_phi1Report else .saem$Gamma2_phi1
  .mat <- nlme::random.effects(.saem)
  .mat2 <- .mat[, .etaTrans, drop = FALSE]
  colnames(.mat2) <- .etaNames
  for (i in seq_along(.eta$name)) {
    .e1 <- .eta$neta1[i]
    .e2 <- .eta$neta2[i]
    .o1 <- .etaTrans[.e1]
    .o2 <- .etaTrans[.e2]
    .ome[.e1, .e2] <- .curOme[.o1, .o2]
    .ome[.e2, .e1] <- .curOme[.o2, .o1]
  }
  env$omega <- .ome
  # Always save the N-row (per-subject) etaMat so FOCEi post-processing
  # gets the correct number of rows, even for mixture models.
  env$.etaMatBase <- .mat2
  if (length(.ui$mixProbs) > 0) {
    .nMix <- length(.ui$mixProbs) + 1
    env$.etaMat <- .mat2[rep(seq_len(nrow(.mat2)), .nMix), , drop = FALSE]
  } else {
    env$.etaMat <- .mat2
  }
  env$etaObf <- data.frame(ID = seq_along(.mat2[, 1]),
                           setNames(as.data.frame(.mat2), .etaNames),
                           OBJI = NA)
  invisible()
}
#' Add Parameter History
#'
#' @param env saem fit environment
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.saemAddParHist <- function(env) {
  .saem <- env$saem
  .ui <- env$ui

  .allThetaNames <- .ui$saemParHistNames

  .m <- .saem$par_hist
  if (ncol(.m) > length(.allThetaNames)) {
    .m <- .m[, seq_along(.allThetaNames)]
  }
  .ph <- data.frame(iter = rep(seq_len(nrow(.m))), as.data.frame(.m),
                    type="Unscaled", check.names=FALSE)
  names(.ph) <- c("iter", .allThetaNames, "type")
  .cls <- class(.ph)
  attr(.cls, "niter") <- env$saemControl$mcmc$niter[1]
  class(.ph) <- .cls
  assign("parHistData", .ph, envir=env)
}
#' Stochastic-approximation (Louis) FIM covariance for SAEM
#'
#' After the estimation iterations, a dedicated covariance phase (`nSaCov`
#' iterations, `src/saem.cpp`) holds the parameters at the converged estimate
#' (gain frozen at 0) and keeps resimulating the individual parameters from their
#' conditional distribution p(phi|y,theta_hat).  The per-iteration Louis
#' observed-information integrand is Monte-Carlo averaged into `saem$HaSa`, giving
#' a converged Fisher information decoupled from the cooling schedule (Monolix's
#' "stochastic approximation" standard errors; Kuhn & Lavielle 2005).  Its inverse
#' is the covariance in (theta, log-Omega-variance, log-sigma2) coordinates; a
#' delta-method Jacobian maps it to the reported scale.
#'
#' The complete-data score currently carries only the diagonal Omega (log-variance)
#' and a single residual variance, so a declared Omega block gets its diagonal
#' variances here and only single-endpoint additive residual SEs are surfaced.
#' @param env saem fit environment
#' @return named full covariance matrix `c(theta, om.<eta>, residual)`, or `NULL`
#' @noRd
.saemSaCov <- function(env) {
  .saemFimToCov(env$saem$HaSa, env)
}
#' Invert a SAEM Fisher Information Matrix into a reported-scale covariance
#'
#' Shared by `covMethod="sa"` (converged FIM `saem$HaSa`) and `covMethod="fim"`
#' (the estimation-phase FIM `saem$Ha`).  Both are the observed information in
#' (theta, log-Omega-variance, log-sigma2) coordinates; this inverts and maps them
#' to the reported scale via a delta-method Jacobian.  The result is required to be
#' positive definite (a noisy/indefinite FIM returns `NULL` so the caller can fall
#' back to the linearized FIM).
#' @param .H Fisher information matrix (nb_param x nb_param)
#' @param env saem fit environment
#' @return named full covariance matrix `c(theta, om.<eta>, residual)`, or `NULL`
#' @noRd
.saemFimToCov <- function(.H, env) {
  .ui <- env$ui
  .saem <- env$saem
  if (is.null(.H) || !is.matrix(.H) || nrow(.H) == 0L ||
        !all(is.finite(.H)) || all(.H == 0)) return(NULL)
  # covariance = inverse of the FIM, in (theta, log-Omega-variance, log-sigma2) coords
  .C <- suppressWarnings(tryCatch(solve(.H), error = function(e) NULL))
  if (is.null(.C) || !all(is.finite(.C))) return(NULL)
  .np <- nrow(.C)
  .tn <- .ui$saemParamsToEstimate[!.ui$saemFixed]
  .nth <- length(.tn)
  if (.nth == 0L || .np < .nth) return(NULL)
  .idf <- .ui$iniDf
  # structural theta block (natural scale; H[1:nth] rows are .tn)
  .ini <- .idf[is.na(.idf$err) & !is.na(.idf$ntheta) & !.idf$fix, "name"]
  if (length(.ui$mixProbs) > 0) .ini <- .ini[!(.ini %in% .ui$mixProbs)]
  .ini <- .ini[.ini %in% .tn]
  .idx <- match(.ini, .tn); .nm <- .ini; .jac <- rep(1, length(.ini))
  # diagonal Omega block: log-variance -> variance, d(var)/d(log var) = var
  .etaN <- tryCatch(.foceiEtaThetaMap(.ui)$etaNames, error = function(e) NULL)
  .omVar <- tryCatch(diag(as.matrix(.saem$Gamma2_phi1)), error = function(e) NULL)
  .nEta <- length(.etaN)
  if (.nEta > 0L && !is.null(.omVar) && length(.omVar) >= .nEta &&
        .np >= .nth + .nEta) {
    .idx <- c(.idx, .nth + seq_len(.nEta))
    .nm <- c(.nm, paste0("om.", .etaN))
    .jac <- c(.jac, .omVar[seq_len(.nEta)])
  }
  # single additive residual: log-sigma2 -> reported SD, d(sd)/d(log sigma2) = 0.5 sd
  .ri <- .idf[!is.na(.idf$err) & !.idf$fix, , drop = FALSE]
  if (nrow(.ri) == 1L && .np == .nth + .nEta + 1L && !grepl("prop|pow", .ri$err)) {
    .ares <- tryCatch(.saem$resMat[1, 1], error = function(e) NA_real_)
    if (is.finite(.ares) && .ares > 0) {
      .idx <- c(.idx, .np); .nm <- c(.nm, .ri$name); .jac <- c(.jac, 0.5 * .ares)
    }
  }
  .cov <- outer(.jac, .jac) * .C[.idx, .idx, drop = FALSE]   # delta method to reported scale
  dimnames(.cov) <- list(.nm, .nm)
  # require a valid (finite, PD) covariance; otherwise let the caller fall back
  if (!all(is.finite(.cov))) return(NULL)
  .ev <- suppressWarnings(tryCatch(eigen(0.5 * (.cov + t(.cov)), symmetric = TRUE,
                                         only.values = TRUE)$values, error = function(e) NA_real_))
  if (any(!is.finite(.ev)) || min(.ev) <= 0) return(NULL)
  .cov
}
#' Splice the linearized-FIM variance block into a fim/sa covariance
#'
#' The analytic (simulation) FIM reliably covers theta + diagonal Omega + additive
#' residuals, but not off-diagonal Omega covariances or proportional/combined residual
#' error (the complete-data Louis correction is unstable when BSV dominates).  For those
#' models this keeps the simulation-based structural-theta block and takes the full
#' variance block (all Omega variances/covariances + residual parameters) from linFim's
#' `calc.COV` (blocB), which handles them correctly via the marginal covariance.  Models
#' the analytic FIM already covers in full are returned unchanged.
#' @param .cov analytic fim/sa covariance (theta + whatever variance params it covers)
#' @param env saem fit environment
#' @return covariance with the linFim variance block spliced in, or `.cov` unchanged
#' @noRd
.saemSpliceLinFimVar <- function(.cov, env) {
  if (!isTRUE(rxode2::rxGetControl(env$ui, "covFull", TRUE))) return(.cov)
  .saem <- env$saem
  attr(.saem, "env") <- env
  .cm <- suppressWarnings(tryCatch(calc.COV(.saem), error = function(e) NULL))
  if (is.null(.cm) || inherits(.cm, "try-error")) return(.cov)
  .vc <- attr(.cm, "varCov")
  if (is.null(.vc) || !is.matrix(.vc) || !all(is.finite(.vc))) return(.cov)
  .vn <- colnames(.vc)
  if (all(.vn %in% rownames(.cov))) return(.cov)     # analytic already covers the variance block
  # keep the simulation structural-theta block; take the whole variance block from linFim
  .rn <- rownames(.cov)
  .th <- .rn[!(.rn %in% .vn) & !grepl("^om\\.|^cov\\.", .rn)]
  .fn <- c(.th, .vn)
  .full <- matrix(0, length(.fn), length(.fn), dimnames = list(.fn, .fn))
  if (length(.th) > 0L) .full[.th, .th] <- .cov[.th, .th, drop = FALSE]
  .full[.vn, .vn] <- .vc
  .full
}

#' Calculate the covariance term
#'
#' @param env saem environment
#' @return nothing, adds $cov to environment if calculated
#' @author Matthew L. Fidler
#' @noRd
.saemCalcCov <- function(env) {
  .ui <- env$ui
  .cm <- rxode2::rxGetControl(.ui, "covMethod", "linFim")
  if (.cm %in% c("sa", "fim")) {
    # Both invert a SAEM observed-information matrix (.saemFimToCov): "sa" uses the
    # converged fixed-theta FIM (saem$HaSa), "fim" the estimation-phase FIM (saem$Ha).
    .H <- if (identical(.cm, "sa")) env$saem$HaSa else env$saem$Ha
    .cov <- NULL
    nlmixrWithTiming("covariance", {
      .cov <- .saemFimToCov(.H, env)
      # off-diagonal Omega / proportional-combined residuals are not reliably in the
      # analytic FIM; splice those from linFim's variance block (blocB).
      if (!is.null(.cov)) .cov <- .saemSpliceLinFimVar(.cov, env)
    })
    if (!is.null(.cov)) {
      # finalization needs a structural-theta cov; stash the full matrix and install
      # it after the fit is built (.saemInstallFullCov).  The control covMethod is reset
      # to its default during finalization, so record the intended label separately.
      .rn <- rownames(.cov)
      .keep <- !grepl("^om\\.|^cov\\.", .rn) & !(.rn %in% .ui$iniDf$name[!is.na(.ui$iniDf$err)])
      env$cov <- .cov[.rn[.keep], .rn[.keep], drop = FALSE]
      assign(".saemFullCov", .cov, envir = env)
      assign(".saemCovMethod", .cm, envir = env)
      env$covMethod <- .cm
      return(invisible())
    }
    message(sprintf("covMethod=\"%s\" could not be computed; using the linearized FIM", .cm))
    rxode2::rxAssignControlValue(.ui, "covMethod", "linFim")
  }
  nlmixrWithTiming("covariance", {
    .ui <- env$ui
    .saem <- env$saem
    attr(.saem, "env") <- env
    .covMethod <- rxode2::rxGetControl(.ui, "covMethod", "linFim")
    .calcCov <- .covMethod == "linFim"
    if (.covMethod == "") {
      .cov <- NULL
      .addCov <- FALSE
    } else {
      .tn <- .ui$saemParamsToEstimate[!.ui$saemFixed]
      .nth <- length(.tn)

      .ini <- .ui$iniDf
      .ini <- .ini[is.na(.ini$err), ]
      .ini <- .ini[!is.na(.ini$ntheta), ]
      .ini <- .ini[!.ini$fix, ]
      .ini <- paste(.ini$name)
      if (length(.ui$mixProbs) > 0) {
        .ini <- .ini[!(.ini %in% .ui$mixProbs)]
      }
      if (.calcCov && .nth == 0) {
        warning("no population parameters in the model, no covariance matrix calculated",
                call.=FALSE)
        .calcCov <- FALSE
        .addCov <- FALSE
        env$cov <- NULL
        .cov <- NULL
        env$covMethod <- "none"
      } else if (.calcCov) {
        .covm <- .saem$Ha[1:.nth, 1:.nth, drop = FALSE]
        ## the FIM linearization (calc.COV) can be ill-conditioned / non-symmetric
        ## (e.g. some delay differential equation models); fail silently and fall
        ## back to the SAEM information matrix rather than aborting the whole fit.
        .covm <- try(calc.COV(.saem), silent = TRUE)
        .doIt <- !inherits(.covm, "try-error")
        if (!.doIt) {
          warning("SAEM covariance by linearization failed; using the SAEM information matrix",
                  call. = FALSE)
        }
        if (.doIt && dim(.covm)[1] != .nth) .doIt <- FALSE
        if (.doIt) {
          # .covm may have NA rows/columns for ill-identified parameters;
          # validate only the well-identified submatrix (.nlmixr2RobustCov()).
          .tmp <- .nlmixr2CholPartial(.covm)
          .addCov <- TRUE
          .sqrtm <- FALSE
          if (inherits(.tmp, "try-error")) {
            .tmp <- .covm
            .tmp <- try(sqrtm(.tmp %*% t(.tmp)), silent = FALSE)
            if (inherits(.tmp, "try-error")) {
              .calcCov <- FALSE
              .covm <- .saem$Ha[1:.nth, 1:.nth, drop = FALSE]
              .tmp <- try(chol(.covm), silent = TRUE)
              .addCov <- TRUE
              .sqrtm <- FALSE
              if (inherits(.tmp, "try-error")) {
                .tmp <- .saem$Ha[1:.nth, 1:.nth, drop = FALSE]
                .tmp <- try(sqrtm(.tmp %*% t(.tmp)), silent = FALSE)
                if (inherits(.tmp, "try-error")) {
                  .addCov <- FALSE
                } else {
                  .sqrtm <- TRUE
                }
              } else {
                .tmp <- .saem$Ha[1:.nth, 1:.nth, drop = FALSE]
              }
            } else {
              .sqrtm <- TRUE
            }
          } else {
            .tmp <- .covm
          }
        } else {
          .tmp <- .saem$Ha[1:.nth, 1:.nth, drop = FALSE]
          .tmp <- try(chol(.tmp), silent = TRUE)
          .calcCov <- FALSE
          .addCov <- TRUE
          .sqrtm <- FALSE
          if (inherits(.tmp, "try-error")) {
            .tmp <- .saem$Ha[1:.nth, 1:.nth, drop = FALSE]
            .tmp <- try(sqrtm(.tmp %*% t(.tmp)), silent = FALSE)
            if (inherits(.tmp, "try-error")) {
              .addCov <- FALSE
            } else {
              .sqrtm <- TRUE
            }
          } else {
            .tmp <- .saem$Ha[1:.nth, 1:.nth, drop = FALSE]
            .calcCov <- FALSE
          }
        }
      } else {
        # non-"linFim" covMethod (0L/"r"/"s"/"r,s"): no calc.COV refinement, use the
        # linearized-FIM Hessian directly (mirrors the calc.COV-failure fallback above).
        .covm <- .saem$Ha[1:.nth, 1:.nth, drop = FALSE]
        .tmp <- try(chol(.covm), silent = TRUE)
        .addCov <- TRUE
        .sqrtm <- FALSE
        if (inherits(.tmp, "try-error")) {
          .tmp <- .saem$Ha[1:.nth, 1:.nth, drop = FALSE]
          .tmp <- try(sqrtm(.tmp %*% t(.tmp)), silent = FALSE)
          if (inherits(.tmp, "try-error")) {
            .addCov <- FALSE
          } else {
            .sqrtm <- TRUE
          }
        } else {
          .tmp <- .saem$Ha[1:.nth, 1:.nth, drop = FALSE]
          .calcCov <- FALSE
        }
      }
      if (.addCov) {
        if (!.calcCov) {
          .cov <- rxode2::rxInv(.tmp)
        } else {
          .cov <- .tmp
        }
        if (!identical(dim(.cov), c(.nth, .nth))) {
          # A degenerate calc.COV (e.g. all parameters unidentified) can make the
          # chol/sqrtm fallback chain silently collapse to the wrong size (unlike
          # calc.COV itself, which is already dimension-checked above); fall back
          # to the linearized-FIM inverse, which is always exactly nth x nth.
          .cov <- tryCatch(rxode2::rxInv(.saem$Ha[1:.nth, 1:.nth, drop = FALSE]),
                            error = function(e) matrix(NA_real_, .nth, .nth))
          .calcCov <- FALSE
        }
        attr(.cov, "dimnames") <- list(.tn, .tn)
        .thCov <- .cov[.ini, .ini, drop = FALSE]           # structural-theta block
        # covFull: assemble the full theta + residual + Omega block-diagonal cov
        # (calc.COV attaches the variance block as "varCov").  The shared output
        # finalization expects a theta-dimensioned cov, so stash the full matrix and
        # install it AFTER the fit is built (.saemInstallFullCov), mirroring focei.
        .vc <- attr(.covm, "varCov")
        .covFull <- isTRUE(rxode2::rxGetControl(.ui, "covFull", TRUE))
        if (.covFull && !is.null(.vc) && is.matrix(.vc) && all(is.finite(.vc))) {
          .vn <- colnames(.vc)
          .fn <- c(.ini, .vn)
          .full <- matrix(0, length(.fn), length(.fn), dimnames = list(.fn, .fn))
          .full[.ini, .ini] <- .thCov
          .full[.vn, .vn] <- .vc
          assign(".saemFullCov", .full, envir = env)
        }
        .cov <- .thCov
      }
    }
    if (.addCov) {
      env$cov <- .cov
      if (.calcCov) {
        env$covMethod <- "linFim"
        if (.addCov & .sqrtm) {
          env$covMethod <- "|linFim|"
          warning("covariance matrix non-positive definite, corrected by sqrtm(linFim %*% linFim)",
                  call.=FALSE)
        }
      } else {
        if (.calcCov) {
          warning("linearization of FIM could not be used to calculate covariance",
                  call.=FALSE)
        }
        if (.addCov & .sqrtm) {
          env$covMethod <- "|fim|"
          warning("covariance matrix non-positive definite, corrected by sqrtm(fim %*% fim)",
                  call.=FALSE)
        } else if (!.addCov) {
          warning("FIM non-positive definite and cannot be used to calculate the covariance",
                  call.=FALSE)
        }
      }
    }

  })
}


#' Get the likelihood name for SAEM
#'
#' @param nnodesGq Number of nodes for Gaussian Quadrature
#' @param nsdGq Number of standard deviations to expand around for
#'   Guassian Quadrature
#' @return Name of the likelihood function
#' @author Matthew L. Fidler
#' @noRd
.saemGetLikName <- function(nnodesGq, nsdGq) {
  if (nnodesGq == 1) {
    paste0("laplace", nsdGq)
  } else {
    paste0("gauss", nnodesGq, "_", nsdGq)
  }
}

#' Calculate the likelihood and the time
#'
#' @param saem saem object
#' @param nnodesGq Number of nodes for Gaussian Quadrature
#' @param nsdGq Number of standard deviations for Gaussian Quadrature
#' @return List of likelihood cauclation time and named objective function value
#' @author Matthew L. Fidler
#' @noRd
.saemCalcLikelihoodTime <- function(saem, nnodesGq, nsdGq, phiM) {
  .likTime <- proc.time()
  .saemObf <- calc.2LL(saem, nnodes.gq = nnodesGq, nsd.gq = nsdGq, phiM)
  .rn <- .saemGetLikName(nnodesGq, nsdGq)
  .likTime <- proc.time() - .likTime
  .likTime <- .likTime["elapsed"]
  list(.likTime, setNames(.saemObf, .rn))
}

#' Calculate the likelihood if requested
#'
#' @param env saem environment
#' @param ...  Other arguments
#' @return Nothing, assigns .saemObf and .likTime in environment
#' @author Matthew L. Fidler
#' @noRd
.saemCalcLikelihood <- function(env, ...) {
  .ui <- env$ui
  env$adjObf <- rxode2::rxGetControl(.ui, "adjObf", TRUE)
  .saem <- env$saem
  .saemCfg <- attr(.saem, "saem.cfg")
  .nphi1 <- .saemCfg$nphi1
  .nphi0 <- .saemCfg$nphi0
  .nphi <- .nphi0 + .nphi1
  .phiM <- matrix(scan(.saemCfg$phiMFile, quiet = TRUE), byrow = TRUE, ncol = .nphi)
  .N <- .saemCfg$N
  dim(.phiM) <- c(.N, .saemCfg$nmc, .saemCfg$niter, .nphi)
  # compresses large object
  env$phiM <- .phiM
  try(unlink(.saemCfg$phiMFile), silent=TRUE)
  .rn <- ""
  .likTime <- 0
  .obf <- rxode2::rxGetControl(.ui, "logLik", FALSE)
  .nnodesGq <- rxode2::rxGetControl(.ui, "nnodesGq", 3)
  .nsdGq <- rxode2::rxGetControl(.ui, "nsd.gq", 1.6)
  if (is.na(.obf)) {
    .saemObf <- NA_real_
  } else if (is.null(.obf)) {
    .tmp <- .saemCalcLikelihoodTime(.saem, .nnodesGq, .nsdGq, env$phiM)
    .likTime <- .tmp[[1]]
    .saemObf <- .tmp[[2]]
  } else if (is(.obf, "logical")) {
    if (is.na(.obf)) {
      .saemObf <- NA_real_
    } else if (.obf) {
      .tmp <- .saemCalcLikelihoodTime(.saem, .nnodesGq, .nsdGq, env$phiM)
      .likTime <- .tmp[[1]]
      .saemObf <- .tmp[[2]]
    } else {
      .saemObf <- NA_real_
    }
  } else if (is(.obf, "numeric")) {
    .saemObf <- .obf
  }
  env$objective <- .saemObf
  env$.likTime  <- .likTime
}

#' Get the calculate cwres residual parameter for saem
#'
#' @param env saem environment
#' @return Calculate resid environment
#' @author Matthew L. Fidler
#' @noRd
.saemGetCalcCwres <- function(env) {
  .ui <- env$ui
  .table <- env$table
  .calcResid <- .table$cwres
  if (is.null(.calcResid)) {
    .calcResid <- .table$saemCWRES
  }
  if (!inherits(.calcResid, "logical")) return(FALSE)
  .calcResid
}

#' Convert the saem options to focei options
#' @param env saem environment that has `$saemControl` assigns focei control to `$control`
#' @return Nothing called for side effects
#' @author Matthew L. Fidler
#' @noRd
.saemControlToFoceiControl <- function(env, assign=TRUE) {
  .saemControl <- env$saemControl
  .ui <- env$ui
  .rxControl <- env$saemControl$rxControl
  # For mixture models the env$.etaMat is the replicated (N*nMix)-row matrix
  # used internally during SAEM.  For the FOCEi post-processing step we need
  # exactly N rows (one per subject).  env$.etaMatBase always holds the
  # N-row version, falling back to env$.etaMat for non-mixture models.
  .etaForFocei <- if (exists(".etaMatBase", envir=env, inherits=FALSE)) {
    env$.etaMatBase
  } else {
    env$.etaMat
  }
  .foceiControl <- foceiControl(maxOuterIterations=0L,
                                maxInnerIterations=0L,
                                covMethod=0L,
                                etaMat=.etaForFocei,
                                sumProd=.saemControl$sumProd,
                                optExpression=.saemControl$optExpression,
                                scaleTo=0,
                                calcTables=.saemControl$calcTables,
                                addProp=.saemControl$addProp,
                                skipCov=.ui$foceiSkipCov,
                                interaction=1L,
                                compress=.saemControl$compress,
                                ci=.saemControl$ci,
                                sigdigTable=.saemControl$sigdigTable,
                                indTolRelax=.saemControl$indTolRelax,
                                rxControl=.rxControl,
                                resetThetaP = 0,
                                resetThetaFinalP = 0,
                                eventSens=.saemControl$eventSens,
                                est = "saem")
  if (exists(".etaMat", envir=env, inherits=FALSE)) {
    rm(list=".etaMat", envir=env)
  }
  if (exists(".etaMatBase", envir=env, inherits=FALSE)) {
    rm(list=".etaMatBase", envir=env)
  }
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @export
#' @rdname nmObjGetFoceiControl
nmObjGetFoceiControl.saem <- function(x, ...) {
  .env <- x[[1]]
  .saemControlToFoceiControl(.env, assign=FALSE)
}

#' Set the extra text for saem
#'
#' @param .env saem environment
#' @param type objective function type
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.setSaemExtra <- function(.env, type) {
  if (inherits(.env, "nlmixr2FitData")) {
    .env <- .env$env
  }
  .ui <- .env$ui
  .txt <- gsub("rxode2 +", "", .ui$modelDesc)
  #.txt <- paste0("(", crayon::italic(ifelse(is.null(.uif$nmodel$lin.solved), ifelse(.uif$predSys, "PRED", "ODE"), "Solved")), "); ")
  .txt <- ""
  if (tolower(type) == "focei") {
    .txt <- paste0(.txt, crayon::silver$italic("OBJF by FOCEi approximation"))
  } else if (tolower(type) == "foce") {
    .txt <- paste0(.txt, crayon::silver$italic("OBJF by FOCE approximation"))
  } else if (tolower(type) == "fo") {
    .txt <- paste0(.txt, crayon::silver$italic("OBJF by FO approximation"))
  } else if (type == "") {
    .txt <- paste0(.txt, crayon::silver$italic("OBJF not calculated"))
  } else {
    .reg <- rex::rex(start, "laplace", capture(.regNum), end)
    .regG <- rex::rex(start, "gauss", capture(.regNum), "_", capture(.regNum), end)
    if (regexpr(.reg, type, perl = TRUE) != -1) {
      .nnode <- 1
      .nsd <- as.numeric(sub(.reg, "\\1", type, perl = TRUE))
    } else if (regexpr(.regG, type, perl = TRUE) != -1) {
      .nnode <- as.numeric(sub(.regG, "\\1", type, perl = TRUE))
      .nsd <- as.numeric(sub(.regG, "\\2", type, perl = TRUE))
    } else {
      stop("unknown error")
    }
    .txt <- paste0(.txt, crayon::silver$italic(sprintf("OBJF by %s", paste0(ifelse(.nnode == 1, "Lapalcian (n.sd=", sprintf("Gaussian Quadrature (n.nodes=%s, n.sd=", .nnode)), .nsd, ")"))))
  }
  .env$extra <- .txt
  invisible()
}

#' Install the stashed full SAEM covariance (theta + residual + Omega) on the fit
#'
#' The shared output finalization expects a theta-dimensioned covariance, so
#' `.saemCalcCov` stashes the full named matrix in `.saemFullCov` and installs the
#' structural-theta block for finalization.  This swaps in the full matrix as
#' `fit$cov` afterward (PD-guarded) and refreshes the parameter table, mirroring
#' focei's `.foceiInstallAnalyticCov`.
#' @param fit saem fit (or its environment)
#' @return nothing, called for side effects
#' @noRd
.saemInstallFullCov <- function(fit) {
  .env <- fit
  if (rxode2::rxIs(fit, "nlmixr2FitData")) .env <- fit$env
  if (!is.environment(.env) || !exists(".saemFullCov", envir = .env, inherits = FALSE)) return(invisible())
  .full <- get(".saemFullCov", envir = .env)
  if (!is.matrix(.full) || !all(is.finite(.full))) return(invisible())
  .full <- 0.5 * (.full + t(.full))                          # exact symmetry (avoids eig_sym warnings)
  .ev <- suppressWarnings(eigen(.full, symmetric = TRUE, only.values = TRUE)$values)
  if (any(diag(.full) <= 0) || !all(is.finite(.ev)) || min(.ev) <= 0) return(invisible())  # keep theta-only
  .env$cov <- .full
  # a valid PD full cov installed: report the intended method (the shared finalization
  # can leave a stale "failed" label even when the SAEM covariance succeeded).  The
  # control covMethod is reset to its default during finalization, so prefer the label
  # recorded by .saemCalcCov (.saemCovMethod) when present.
  .m <- if (exists(".saemCovMethod", envir = .env, inherits = FALSE)) {
    get(".saemCovMethod", envir = .env)
  } else {
    tryCatch(rxode2::rxGetControl(.env$ui, "covMethod", "linFim"), error = function(e) "linFim")
  }
  .env$covMethod <- if (.m %in% c("sa", "fim")) .m else "linFim"
  # surface the residual (error-model theta) SEs in the parameter table from the full
  # cov -- these are theta rows with a missing SE (Omega variances are reported as BSV,
  # with their SEs available in $cov).
  if (exists("parFixedDf", envir = .env, inherits = FALSE)) {
    .pf <- .env$parFixedDf
    .se <- sqrt(diag(.full))
    .ci <- tryCatch(as.numeric(rxode2::rxGetControl(.env$ui, "ci", 0.95)), error = function(e) 0.95)
    .qn <- stats::qnorm(1 - (1 - .ci) / 2)
    for (.n in rownames(.pf)) {
      if (.n %in% names(.se) && "SE" %in% names(.pf) &&
            (is.na(.pf[.n, "SE"]) || !is.finite(.pf[.n, "SE"]) || .pf[.n, "SE"] < 1e-100)) {
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
    .env$parFixedDf <- .pf
  }
  invisible()
}

#' Fit the saem family of models
#'
#' @param env Environment from nlmixr2Est
#' @param ... Other arguments
#' @return fit environment with $saem, $saemControl, $dataSav, $origData, $ui
#' @author Matthew L. Fidler
#' @noRd
.saemFamilyFit <- function(env, ...) {
  .ui <- env$ui
  .control <- .ui$control
  .data <- env$data
  .ret <- new.env(parent=emptyenv())
  .ret$table <- env$table
  nlmixrWithTiming("setup", {
    .foceiPreProcessData(.data, .ret, .ui, .control$rxControl)
    .tv <- .nlmixrTimeVaryingCovariates(.ret$dataSav, .ui, .control$rxControl)
  })

  .ret$saem <- .saemFitModel(.ui, .ret$dataSav, timeVaryingCovariates=.tv)
  # Re-stage the mu-ref time-varying split for the post-processing: the theta
  # table and parameter history are named from saemParamsToEstimate/
  # saemParHistNames, which only put a time-varying covariate in the correct
  # (Plambda) order while the split is staged (.saemFitModel stages it only for
  # the fit itself, then unstages on exit).
  .nlmixrSetMuRefTimeVarying(.ui, .tv)
  on.exit(.nlmixrRmMuRefTimeVarying(.ui), add = TRUE)
  .ret$ui <- .ui
  .saemCalcCov(.ret)
  .ret <- nlmixrWithTiming("postprocess", {
    if (!is.null(.ret$saem$tolFactor)) {
      .ret$tolFactor <- .ret$saem$tolFactor
    }
    .ret$control <- .control
    nmObjHandleControlObject(.ret$control, .ret)
    .getSaemTheta(.ret)
    .getSaemOmega(.ret)
    # Must run against the un-pooled omega, before .saemMixFix() pools split
    # ETAs, or ui$theta silently falls back to ini() values for every param.
    .nlmixr2FitUpdateParams(.ret)
    # Builds mixList/mixNum/mixIcov; must run before nlmixr2CreateOutputFromUi.
    .saemMixFix(.ret, .ui)
    .ui <- .ret$ui
    .saemAddParHist(.ret)
    .saemCalcLikelihood(.ret)
    if (is.environment(.ui) && exists("control", envir=.ui, inherits=FALSE)) {
      rm(list="control", envir=.ui)
    }
    .ret$theta <- .ui$saemThetaDataFrame
    .ret$model <- .ui$saemModelPred
    .ret$message <- "" # no message for now
    .ret$est <- "saem"
    .saemControlToFoceiControl(.ret)
    .ret <- nlmixr2CreateOutputFromUi(.ret$ui, data=.ret$origData, control=.ret$control, table=.ret$table, env=.ret, est="saem")
    # covFull/sa: swap in the stashed full theta+residual+Omega covariance now that
    # the theta-dimensioned fit table has been built.
    .saemInstallFullCov(.ret)
    # The shared output finalization can leave a stale "failed" covMethod label even
    # when the SAEM covariance actually succeeded (a valid finite PD $cov exists);
    # restore the intended method in that case.
    .rEnv <- if (rxode2::rxIs(.ret, "nlmixr2FitData")) .ret$env else .ret
    if (is.environment(.rEnv) && identical(.rEnv$covMethod, "failed") &&
          is.matrix(.rEnv$cov) && all(is.finite(.rEnv$cov)) && all(diag(.rEnv$cov) > 0)) {
      # the control covMethod is reset to its default during finalization, so prefer the
      # label recorded by .saemCalcCov (.saemCovMethod) when present.
      .cm <- if (exists(".saemCovMethod", envir = .rEnv, inherits = FALSE)) {
        get(".saemCovMethod", envir = .rEnv)
      } else {
        tryCatch(rxode2::rxGetControl(.ui, "covMethod", "linFim"), error = function(e) "linFim")
      }
      .rEnv$covMethod <- if (.cm %in% c("linFim", "fim", "sa")) .cm else "linFim"
    }
    # For mixture models: post-correct me/mn/mu in the assembled fit table
    # (mirrors the .mixFixTable call in .foceiFamilyReturn for FOCEi fits)
    if (inherits(.ret, "nlmixr2FitData") && length(.ui$mixProbs) > 0L) {
      .retEnv <- attr(class(.ret), ".foceiEnv")
      if (is.null(.retEnv)) .retEnv <- .ret$env
      .ret <- .mixFixTable(.ret, .retEnv, .ui)
    }
    .setSaemExtra(.ret, "FOCEi")
    .ret
  })
  .env <- .ret$env
  .env$method <- "SAEM "
  .ret
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.saem <- function(env, ...) {
  .ui <- env$ui
  # saem supports a general log-likelihood endpoint (ll() ~ expr) the same way
  # saemix does (the model returns the per-obs loglik; the RWM kernels use -ll as
  # the observation loss); only require normality for the ordinary case.
  if (!.fsaemGeneralLik(.ui)) {
    rxode2::assertRxUiTransformNormal(.ui, " for the estimation routine 'saem'", .var.name=.ui$modelName)
  }
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'saem'",
                             .var.name=.ui$modelName)
  rxode2::assertRxUiMixedOnly(.ui, " for the estimation routine 'saem'", .var.name=.ui$modelName)
  rxode2::warnRxBounded(.ui, " which are ignored in 'saem'", .var.name=.ui$modelName)
  if (length(.ui$mixProbs) > 0) {
    message("mixture SAEM computation scales with the number of sub-populations")
  }
  .saemFamilyControl(env, ...)
  on.exit({
    if (is.environment(.ui) && exists("control", envir=.ui, inherits=FALSE)) {
      rm("control", envir=.ui)
    }
  }, add=TRUE)
  .saemFamilyFit(env,  ...)
}
attr(nlmixr2Est.saem, "covPresent") <- TRUE
attr(nlmixr2Est.saem, "unbounded") <- TRUE
attr(nlmixr2Est.saem, "mu") <- TRUE
attr(nlmixr2Est.saem, "iov") <- TRUE


#' @rdname nmObjGet
#' @export
nmObjGet.saemDopredIpred <- function(x, ...) {
  .obj <- x[[1]]
  .env <- .obj$env
  if (exists("saem0", envir=.env)) {
    .saem <- .env$saem
    .saemCfg <- attr(.saem, "saem.cfg")
    .dopred <- attr(.saem, "dopred")
    .dopred(.saem$mpost_phi, .saemCfg$evt, .saemCfg$opt)
  } else {
    NULL
  }
}

#' @rdname nmObjGet
#' @export
nmObjGet.saemDopredPred <- function(x, ...) {
  .obj <- x[[1]]
  .env <- .obj$env
  if (exists("saem0", envir=.env)) {
    .saem <- .env$saem
    .saemCfg <- attr(.saem, "saem.cfg")
    .dopred <- attr(.saem, "dopred")
    .dopred(.saem$mprior_phi, .saemCfg$evt, .saemCfg$opt)
  } else {
    NULL
  }
}
#attr(nmObjGet.saemDopredIpred, "desc") <- "Get ipred from low level saem"
