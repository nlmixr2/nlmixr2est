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
  .muRefCovariateDataFrame <- ui$muRefCovariateDataFrame
  if (length(timeVaryingCovariates) > 0) {
    # Drop time-varying covariates
    # First get the time varying covariates
    .w <- which(.muRefCovariateDataFrame$covariate %in% timeVaryingCovariates)
    # next find out the theta for the phi expression
    .covPar <- .muRefCovariateDataFrame[.w, "theta"]
    .w2 <- which(ui$muRefCurEval$parameter %in% .covPar)
    if (length(.w2) > 0) {
      # see if the expression is on a log scale
      .w3 <- which("exp" == ui$muRefCurEval$curEval[.w2])
      if (length(.w3) > 0) {
        .w2 <- .w2[.w3]
        .texp <- ui$muRefCurEval$parameter[.w2]
        # now get parameters
        .pars <- .muRefCovariateDataFrame$covariateParameter[.muRefCovariateDataFrame$theta %in% .texp]
        ## warning(paste0("log-scale mu referenced time varying covariates (",
        ##                paste(.pars, collapse=", "),
        ##                ") may have better results on no log-transformed scale (https://github.com/nlmixr2/nlmixr2est/issues/348), check results for plausibility"),
        ##         call.=FALSE)
      }

    }
    .muRefCovariateDataFrame <- .muRefCovariateDataFrame[!(.muRefCovariateDataFrame$covariate %in% timeVaryingCovariates), ]
  }
  assign("muRefFinal", .muRefCovariateDataFrame, ui)
  assign("timeVaryingCovariates", timeVaryingCovariates, ui)
  on.exit({
    if (is.environment(ui) && exists("muRefFinal", envir=ui, inherits=FALSE)) {
      rm(list="muRefFinal", envir=ui)
    }
    if (is.environment(ui) && exists("timeVaryingCovariates", envir=ui, inherits=FALSE)) {
      rm(list="timeVaryingCovariates", envir=ui)
    }
  })
  .model <- ui$saemModelList
  .inits <- ui$saemInit
  .rxControl <- rxode2::rxGetControl(ui, "rxControl", rxode2::rxControl())
  .ue <- .uninformativeEtas(ui,
                            handleUninformativeEtas=rxode2::rxGetControl(ui, "handleUninformativeEtas", TRUE),
                            data=data,
                            attr(.model$saem_mod, "rx"),
                            rxControl=.rxControl)
  .cfg <- nlmixrWithTiming("configure", {
    .cfg <- .configsaem(model=.model,
                        data=data,
                        inits=.inits,
                        mcmc=rxode2::rxGetControl(ui, "mcmc",
                                                  list(niter = c(200, 300),
                                                       nmc = 3, nu = c(2, 2, 2))),
                        rxControl=.rxControl,
                        distribution="normal",
                        fixedOmega=ui$saemModelOmegaFixed,
                        fixedOmegaValues=ui$saemModelOmegaFixedValues,
                        parHistThetaKeep=ui$saemParHistThetaKeep,
                        parHistOmegaKeep=ui$saemParHistOmegaKeep,
                        seed=rxode2::rxGetControl(ui, "seed", 99),
                        DEBUG=rxode2::rxGetControl(ui, "DEBUG", 0),
                        tol=rxode2::rxGetControl(ui, "tol", 1e-6),
                        itmax=rxode2::rxGetControl(ui, "itmax", 30),
                        type=rxode2::rxGetControl(ui, "type", "nelder-mead"),
                        lambdaRange=rxode2::rxGetControl(ui, "lambdaRange", 3),
                        powRange=rxode2::rxGetControl(ui, "powRange", 10),
                        odeRecalcFactor=rxode2::rxGetControl(ui, "odeRecalcFactor", 10^0.5),
                        maxOdeRecalc=rxode2::rxGetControl(ui, "maxOdeRecalc", 10^0.5),
                        indTolRelax=rxode2::rxGetControl(ui, "indTolRelax", TRUE),
                        nres=ui$saemModNumEst,
                        perSa=rxode2::rxGetControl(ui, "perSa", 0.75),
                        perNoCor=rxode2::rxGetControl(ui, "perNoCor", 0.75),
                        perFixOmega=rxode2::rxGetControl(ui, "perFixOmega", 0.1),
                        perFixResid=rxode2::rxGetControl(ui, "perFixResid", 0.1),
                        resFixed=ui$saemResFixed,
                        ue=.ue,
                        mixProb=ui$saemMixProb,
                        mixProbMethod=rxode2::rxGetControl(ui, "mixProbMethod", "regularized"),
                        mixProbStepExp=rxode2::rxGetControl(ui, "mixProbStepExp", 1),
                        mixProbPriorN=rxode2::rxGetControl(ui, "mixProbPriorN", 20),
                        mixSampleMethod=rxode2::rxGetControl(ui, "mixSampleMethod", "parallel"),
                        omegaShare=ui$saemOmegaShare,
                        omegaShareSubpop=ui$saemOmegaShareSubpop)
    .cfg$cres <- ui$saemCres
    .cfg$yj <- ui$saemYj
    .cfg$lres <- ui$saemLres
    .cfg$low <- ui$saemLow
    .cfg$hi <- ui$saemHi
    .cfg$propT <- ui$saemPropT
    .cfg$addProp <- ui$saemAddProp
    .cfg$resValue <- ui$saemResValue
    # Iteration-print flows through the shared scaleApplyIterPrintControl
    # (src/scale.h); U auto-skips since Plambda has no optimizer scaling,
    # leaving # and X. xform (xPar/probitIdx/bounds) drives the X row's
    # back-transform via scaleAttachXform, same as every other estimator.
    .cfg$parHistNames <- as.character(ui$saemParHistNames)
    .cfg$xform        <- .iterPrintXParFromUi(ui, .cfg$parHistNames)
    .cfg$iterPrintControl <- rxode2::rxGetControl(ui, "iterPrintControl",
                                                  iterPrintControl())
    .saemCheckCfg(.cfg)
    .cfg
  })
  nlmixrWithTiming("saem", {
    .model$saem_mod(.cfg)
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
#' Online batch-means stochastic-approximation covariance for SAEM
#'
#' The averaged SAEM parameter iterates are asymptotically normal, so the
#' asymptotic covariance is estimated directly from the recorded cooling-phase
#' iterate trajectory (`saem$par_hist`) by the online batch-means estimator of
#' Jiang, Roy, Balasubramanian, Davis, Drusvyatskiy & Na (2025).  It is
#' derivative-free and inherently full -- theta + Omega (diagonal) variance +
#' residual -- since every one of those is an iterated parameter recorded on its
#' reported scale.  Off-diagonal Omega covariances are not in `par_hist`, so a
#' declared block gets only its diagonal variances here.
#' @param env saem fit environment
#' @return named full covariance matrix `c(theta, om.<eta>, residual)`, or `NULL`
#' @noRd
.saemSaCov <- function(env) {
  .ui <- env$ui
  .saem <- env$saem
  .ph <- .saem$par_hist
  if (is.null(.ph) || !is.matrix(.ph)) return(NULL)
  .nm <- as.character(.ui$saemParHistNames)
  if (ncol(.ph) < length(.nm)) return(NULL)
  .ph <- .ph[, seq_along(.nm), drop = FALSE]
  colnames(.ph) <- .nm
  # drop mixture-probability columns -- not part of the parameter covariance
  .mix <- .ui$mixProbs
  if (length(.mix) > 0L) .ph <- .ph[, !(.nm %in% .mix), drop = FALSE]
  # cooling / smoothing phase = the second SAEM phase (nEm iterations after nBurn)
  .niter <- tryCatch(env$saemControl$mcmc$niter, error = function(e) NULL)  # c(nBurn, nEm)
  .N <- nrow(.ph)
  .nBurn <- if (length(.niter) >= 1L && is.finite(.niter[1])) as.integer(.niter[1]) else as.integer(floor(.N / 2))
  .cool <- if (.N > .nBurn && .nBurn >= 0L) seq.int(.nBurn + 1L, .N) else seq_len(.N)
  X <- .ph[.cool, , drop = FALSE]
  .n <- nrow(X)
  if (.n < 20L || ncol(X) == 0L) return(NULL)           # too few cooling iterations to be meaningful
  Xc <- sweep(X, 2, colMeans(X))                        # centered iterates
  # increasing block schedule a_m = floor(C m^beta) (Jiang et al. 2025, Remark 1)
  .C <- suppressWarnings(as.numeric(rxode2::rxGetControl(.ui, "covSaBlock", 1)))
  if (!is.finite(.C) || .C < 1) .C <- 1
  .beta <- 2
  .breaks <- unique(pmin(.n, floor(.C * (seq_len(.n))^.beta)))
  .breaks <- c(0L, .breaks[.breaks >= 1L])
  if (utils::tail(.breaks, 1L) != .n) .breaks <- c(.breaks, .n)
  .breaks <- unique(.breaks)
  .Sig <- matrix(0, ncol(X), ncol(X))
  for (.i in seq_len(length(.breaks) - 1L)) {
    .idx <- seq.int(.breaks[.i] + 1L, .breaks[.i + 1L])
    .Sm <- colSums(Xc[.idx, , drop = FALSE])
    .Sig <- .Sig + tcrossprod(.Sm)
  }
  # Sigmahat = sum_m S_m S_m' / n estimates the limiting cov of sqrt(n)*(xbar - x*);
  # the covariance of the averaged estimate xbar is Sigmahat / n.
  .cov <- .Sig / (.n * .n)
  # name: theta as-is; V(<eta>) -> om.<eta> (diagonal Omega variance); residual as-is
  .rn <- colnames(X)
  .isV <- grepl("^V\\(", .rn)
  .rn[.isV] <- paste0("om.", sub("^V\\((.*)\\)$", "\\1", .rn[.isV]))
  dimnames(.cov) <- list(.rn, .rn)
  .cov
}

#' Calculate the covariance term
#'
#' @param env saem environment
#' @return nothing, adds $cov to environment if calculated
#' @author Matthew L. Fidler
#' @noRd
.saemCalcCov <- function(env) {
  .ui <- env$ui
  if (identical(rxode2::rxGetControl(.ui, "covMethod", "linFim"), "sa")) {
    # The online batch-means SA covariance (.saemSaCov) is implemented but requires a
    # dedicated slow-decaying-gain phase to match the estimate covariance -- SAEM's
    # default gain schedule collapses the cooling trajectory and underestimates it.
    # Until that phase is wired up, fall back to the linearized FIM.
    message("covMethod=\"sa\" (online batch-means) is not yet available; using the linearized FIM")
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
        .covm <- try(calc.COV(.saem))
        .doIt <- !inherits(.covm, "try-error")
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
          .tmp <- try(chol(.covm), silent = TRUE)
          .tmp <- .saem$Ha[1:.nth, 1:.nth]
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
  .ev <- suppressWarnings(eigen(.full, symmetric = TRUE, only.values = TRUE)$values)
  if (any(diag(.full) <= 0) || !all(is.finite(.ev)) || min(.ev) <= 0) return(invisible())  # keep theta-only
  .env$cov <- .full
  # a valid PD full cov installed: report the intended method (the shared finalization
  # can leave a stale "failed" label even when the SAEM covariance succeeded)
  .m <- tryCatch(rxode2::rxGetControl(.env$ui, "covMethod", "linFim"), error = function(e) "linFim")
  .env$covMethod <- if (identical(.m, "sa")) "sa" else "linFim"
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
  rxode2::assertRxUiTransformNormal(.ui, " for the estimation routine 'saem'", .var.name=.ui$modelName)
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
