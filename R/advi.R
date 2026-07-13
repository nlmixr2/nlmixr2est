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

#' Run the ADVI optimization: prep, inner setup, initialize the variational +
#' population state, and drive the whole optimization (the adaptive step-size
#' search + the main loop) in one C++ call (adviOptimize_).
#' @param ui bounded-transformed rxode2 ui
#' @param data estimation data
#' @param control adviControl
#' @param resume optional list from a previous fit's `$adviState` for warm resume
#' @return the raw ADVI result list (variational params, estimates, elbo, parHist)
#' @noRd
.adviOptimize <- function(ui, data, control, resume = NULL) {
  .prep <- .adviDataPrep(ui, data)
  N <- .prep$N; neta <- .prep$neta
  ## a resumed run keeps its original family; otherwise use the control's
  .fr <- if (!is.null(resume) && !is.null(resume$family))
    identical(resume$family, "fullRank") else identical(control$adviFamily, "fullRank")

  ## the FOCEi inner setup starts at the variational means (resumed or 0); the
  ## optimization state itself is initialized/resumed inside adviOptimize_
  .etaMat0 <- if (is.null(resume)) matrix(0, N, neta) else resume$mu
  .setup <- .adviInnerSetup(ui, data, .etaMat0, control)
  on.exit(.adviInnerFree(), add = TRUE)

  ## iteration printing: the shared scale.h table (like saem/vae).  Rows are
  ## always captured (-> standard parHistData); iterPrintControl$every gates the
  ## console output.  The step-size search runs join the same table as labeled
  ## "srch <eta>" phases; the main run is the "SGA" phase.
  .ipNames <- c(.prep$thetaRealNames, paste0("o(", .prep$etaNames, ")"))
  .ipXform <- .iterPrintXParFromUi(rxode2::rxUiDecompress(ui), .ipNames)

  ## thread count for the parallel per-subject ELBO core (same knob as the inner
  ## eval driver: rxControl$cores, falling back to the rxode2 thread pool).  Kept
  ## bit-for-bit invariant to the thread count by a serial id-ordered reduction.
  .cores <- tryCatch({
    .c <- control$rxControl$cores
    if (is.null(.c) || is.na(.c) || .c < 1L) as.integer(rxode2::getRxThreads()) else as.integer(.c)
  }, error = function(e) 1L)

  ## everything else -- state init/resume, the mu-ref and full-Bayes phi maps,
  ## the adaptEta search, the main loop, and the derived result fields -- runs
  ## in one C++ call (a resumed seed/etaScale is picked up from `resume`)
  .res <- adviOptimize_(list(
    pointEstimate = isTRUE(control$pointEstimate), fr = as.integer(.fr),
    N = as.integer(N),
    theta = as.numeric(.prep$theta), omega = as.numeric(.prep$omega),
    muRefThetaIdx = as.integer(.prep$muRefThetaIdx),
    thetaFix = as.logical(.prep$thetaFix), omegaFix = as.logical(.prep$omegaFix),
    iters = as.integer(control$iters), seed = as.numeric(control$seed),
    tau = as.numeric(control$tau), alpha = as.numeric(control$alpha),
    nMc = as.integer(control$nMc), cores = .cores,
    adaptEta = isTRUE(control$adaptEta),
    etaCandidates = as.numeric(control$etaCandidates),
    nAdapt = as.integer(min(control$iters, 75L)),
    parNames = .ipNames, iterPrintControl = control$iterPrintControl,
    xform = .ipXform, resume = resume))
  .res$family <- control$adviFamily
  .res$prep <- .prep
  .res$etaNames <- .prep$etaNames
  .res$thetaNames <- names(.prep$th)
  .res$model <- .setup$model
  class(.res) <- "nlmixr2advi"
  .res
}

#' Assemble the standard nlmixr2FitData from an ADVI result: seed the ui iniDf
#' with the ADVI estimates (population thetas + between-subject omega diagonal),
#' supply the variational posterior means as the FOCEi inner EBE start (etaMat),
#' and run the eval-only FOCEi finalize (maxOuterIterations=0) which reuses
#' inner.cpp for the objective, EBEs, residual tables, and the covariance step.
#' No outer optimizer is run; the ADVI estimates are final.  Mirrors .vaeToFit /
#' .rpemBuildFit.
#' @noRd
.adviToFit <- function(env, res) {
  .ui <- env$ui
  .control <- env$adviControl
  .prep <- res$prep
  .rxControl <- .control$rxControl

  .ret <- new.env(parent = emptyenv())
  .ret$table <- env$table
  .foceiPreProcessData(env$data, .ret, .ui, .rxControl)

  ## seed the ui iniDf with the ADVI estimates so the eval reports them
  .uiD <- rxode2::rxUiDecompress(.ui)
  .idf <- .uiD$iniDf
  .thRow <- !is.na(.idf$ntheta)
  .idf$est[.thRow] <- res$theta[.idf$ntheta[.thRow]]
  .popOm <- stats::setNames(res$popOmega, .prep$etaNames)
  .etaRow <- !is.na(.idf$neta1) & .idf$neta1 == .idf$neta2
  .idf$est[.etaRow] <- .popOm[.idf$name[.etaRow]]
  assign("iniDf", .idf, envir = .uiD)
  .ui2 <- rxode2::rxUiCompress(.uiD)

  ## variational posterior means as the FOCEi inner EBE start [nsub, neta]
  .eb <- res$mu
  colnames(.eb) <- .prep$etaNames
  .ret$.etaMat <- .eb
  .ret$.etaMatBase <- .eb
  .ret$etaObf <- data.frame(ID = seq_len(nrow(.eb)),
                            stats::setNames(as.data.frame(.eb), .prep$etaNames),
                            OBJI = NA)
  .om <- diag(res$popOmega, .prep$neta)
  dimnames(.om) <- list(.prep$etaNames, .prep$etaNames)
  .ret$omega <- .om
  .ret$ui <- .ui2
  .ret$fullTheta <- stats::setNames(res$theta, names(.prep$th))

  ## covMethod="advi": for full-Bayes the SEs come from the population variational
  ## covariance (installed below, so skip the FOCEi cov step); for point-estimate
  ## there is no population variational block, so fall back to the FOCEi "r,s".
  .covM <- if (identical(.control$covMethod, "advi"))
    (if (isTRUE(res$pointEstimate)) "r,s" else "") else .control$covMethod
  .lik <- .control$likelihood
  .interaction <- if (.lik %in% c("foce", "focep")) 0L else 1L
  .foce <- if (identical(.lik, "focep")) "foce+" else "nonmem"
  .fc <- foceiControl(rxControl = .rxControl, maxOuterIterations = 0L,
                      maxInnerIterations = 0L, covMethod = .covM, etaMat = .eb,
                      scaleTo = 0, interaction = .interaction, foce = .foce,
                      sumProd = .control$sumProd, optExpression = .control$optExpression,
                      literalFix = .control$literalFix, literalFixRes = .control$literalFixRes,
                      addProp = .control$addProp, calcTables = .control$calcTables,
                      compress = .control$compress, ci = .control$ci,
                      sigdigTable = .control$sigdigTable, stickyRecalcN = .control$stickyRecalcN,
                      maxOdeRecalc = .control$maxOdeRecalc, odeRecalcFactor = .control$odeRecalcFactor,
                      indTolRelax = .control$indTolRelax, eventSens = .control$eventSens,
                      fast = FALSE, print = 0L)
  .ret$control <- .fc
  .ret$method <- "advi"
  .ret$extra <- ""
  .ret$est <- "advi"
  .ret$ofvType <- "advi"
  .ret$adjObf <- .control$adjObf
  ## the ADVI optimization walk (standard parHistData -> $parHist accessor)
  if (!is.null(res$parHistData)) .ret$parHistData <- res$parHistData
  nmObjHandleControlObject(.control, .ret)   # store adviControl for nmObjGetControl.advi
  ## reuse the models compiled for the ADVI loop (inner/EBE/pred + thetaSens):
  ## with $model present the eval-only finalize skips its own symengine rebuild
  ## (the finalize reads only the foce-prefix columns of the inner model, so the
  ## interaction-model column layout is compatible)
  if (!is.null(res$model)) {
    .ret$model <- res$model
  } else {
    .ret$foceiModel <- .ui2$focei
  }
  .fit <- nlmixr2CreateOutputFromUi(.ret$ui, data = .ret$origData, control = .fc,
                                    table = .ret$table, env = .ret, est = "advi")
  ## ADVI artifacts + warm-resume state on the fit env
  .e <- .fit$env
  .e$adviElbo <- res$elbo
  .st <- list(mu = res$mu, theta = res$theta, logPopOmega = res$logPopOmega,
              it0 = res$it0, sMu = res$sMu, sScale = res$sScale, sTheta = res$sTheta,
              sLpo = res$sLpo, seed = res$seed, etaScale = res$etaScale,
              family = res$family, pointEstimate = res$pointEstimate)
  if (isTRUE(res$pointEstimate)) {
    if (identical(res$family, "fullRank")) .st$Lpack <- res$scale else .st$omega <- res$scale
  } else {
    ## full-Bayes: per-subject scale is generic; also persist the population block
    .st$scale <- res$scale; .st$mPop <- res$mPop; .st$Lpop <- res$Lpop
    .st$smPop <- res$smPop; .st$sLpop <- res$sLpop
    ## population variational covariance -> named phi-space cov on the fit env
    .cov <- res$adviCov
    .thNm <- res$prep$thetaRealNames[res$phiThetaIdx[res$phiThetaIdx >= 0] + 1L]
    .nm <- c(.thNm, paste0("omega.", res$etaNames[res$phiOmIdx[res$phiOmIdx >= 0] + 1L]))
    if (nrow(.cov) == length(.nm)) dimnames(.cov) <- list(.nm, .nm)
    .e$adviCov <- .cov
    ## install the population variational covariance as the fit's SE source (the
    ## theta block maps directly to parFixedDf's population/residual parameters)
    .adviInstallVarCov(.fit, res)
  }
  .e$adviState <- .st
  .fit
}

#' Install the ADVI population variational covariance (Lpop Lpop^T) as the fit's
#' covariance + parFixedDf SEs -- the natural full-Bayes uncertainty.  The theta
#' block of the phi-space covariance maps directly to the population / residual
#' parameters (by name); the log-variance block is retained on $env$adviCov.
#' Reuses the rpem Fisher-cov installer.
#' @noRd
.adviInstallVarCov <- function(fit, res) {
  .thComp <- which(res$phiThetaIdx >= 0)
  if (length(.thComp) == 0L) return(invisible())
  .thNames <- res$prep$thetaRealNames[res$phiThetaIdx[.thComp] + 1L]
  .thetaCov <- res$adviCov[.thComp, .thComp, drop = FALSE]
  dimnames(.thetaCov) <- list(.thNames, .thNames)
  .rpemInstallFisherCov(fit, .thetaCov)
  .env <- if (rxode2::rxIs(fit, "nlmixr2FitData")) fit$env else fit
  if (is.environment(.env)) .env$covMethod <- "advi"
  invisible()
}

#' Fit an ADVI model: set up the inner/outer problems and run the C++ loop.
#' @param env estimation environment (holds ui, data, adviControl)
#' @noRd
.adviFitModel <- function(env) {
  .ui <- env$ui
  .control <- env$adviControl
  ## warm resume: accept a prior advi fit or its adviState
  .resume <- .control$resume
  if (!is.null(.resume)) {
    if (rxode2::rxIs(.resume, "nlmixr2FitData")) .resume <- .resume$env$adviState
    else if (is.environment(.resume) && exists("adviState", .resume)) .resume <- .resume$adviState
    if (!is.list(.resume) || is.null(.resume$it0))
      stop("est=\"advi\" 'resume' must be a prior advi fit or its $env$adviState", call. = FALSE)
  }
  .res <- .adviOptimize(.ui, env$data, .control, resume = .resume)
  if (isTRUE(.control$returnAdvi)) return(.res)
  .adviToFit(env, .res)
}
