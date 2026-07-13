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

#' Adaptively choose the step-size scale `eta` (paper Sec 2.6): run a short loop
#' (from the initial state) for each candidate and pick the one with the best
#' late-iteration mean ELBO.  Deterministic (same seed/init), so it does not
#' break the prefix/resume reproducibility of the main run.
#' @noRd
.adviAdaptEta <- function(mu0, omega0, theta0, logPopOmega0, muRefIdx, thetaFix,
                          omegaFix, control, seed) {
  .cands <- control$etaCandidates
  if (length(.cands) == 1L) return(.cands)
  N <- nrow(mu0); neta <- ncol(mu0); ntheta <- length(theta0)
  .nAdapt <- as.integer(min(control$iters, 30L))
  .z <- function() matrix(0, N, neta)
  .best <- -Inf; .bestEta <- .cands[1]
  for (.e in .cands) {
    .r <- adviLoop_(mu0, omega0, theta0, logPopOmega0, muRefIdx, thetaFix, omegaFix,
                    .nAdapt, as.numeric(seed), .e, as.numeric(control$tau),
                    as.numeric(control$alpha), as.integer(control$nMc), 0L,
                    .z(), .z(), numeric(ntheta), numeric(neta))
    .el <- .r$elbo
    .score <- if (all(is.finite(.el)))
      mean(.el[max(1L, length(.el) - .nAdapt %/% 3L):length(.el)]) else -Inf
    if (.score > .best) { .best <- .score; .bestEta <- .e }
  }
  .bestEta
}

#' Run the ADVI optimization: prep, inner setup, initialize the variational +
#' population state, and drive the 100%-C++ loop (adviLoop_).
#' @param ui bounded-transformed rxode2 ui
#' @param data estimation data
#' @param control adviControl
#' @param resume optional list from a previous fit's `$adviState` for warm resume
#' @return the raw ADVI result list (variational params, estimates, elbo, parHist)
#' @noRd
.adviOptimize <- function(ui, data, control, resume = NULL) {
  .prep <- .adviDataPrep(ui, data)
  N <- .prep$N; neta <- .prep$neta; ntheta <- .prep$ntheta

  ## initial variational + population state (or resume from a prior fit)
  if (is.null(resume)) {
    .logPopOmega0 <- log(.prep$omega)
    .mu0 <- matrix(0, N, neta)
    ## start q at the prior scale: log-sd = 0.5 log(popOmega)
    .omega0 <- matrix(rep(0.5 * .logPopOmega0, each = N), N, neta)
    .theta0 <- .prep$theta
    .it0 <- 0L
    .sMu <- matrix(0, N, neta); .sOmega <- matrix(0, N, neta)
    .sTheta <- numeric(ntheta); .sLpo <- numeric(neta)
  } else {
    .mu0 <- resume$mu; .omega0 <- resume$omega
    .theta0 <- resume$theta; .logPopOmega0 <- resume$logPopOmega
    .it0 <- as.integer(resume$it0)
    .sMu <- resume$sMu; .sOmega <- resume$sOmega
    .sTheta <- resume$sTheta; .sLpo <- resume$sLpo
  }

  .adviInnerSetup(ui, data, .mu0, control)
  on.exit(.adviInnerFree(), add = TRUE)

  ## the counter-based RNG is keyed by the global iteration index, so resuming
  ## with the original seed continues the exact same stream (prefix property).
  .seed <- if (!is.null(resume) && !is.null(resume$seed)) resume$seed else control$seed
  ## step-size scale: reuse the resumed run's, else adaptively search, else fixed.
  .muRefIdx <- as.integer(.prep$muRefThetaIdx)
  .thFix <- as.logical(.prep$thetaFix); .omFix <- as.logical(.prep$omegaFix)
  .etaScale <- if (!is.null(resume) && !is.null(resume$etaScale)) resume$etaScale
    else if (isTRUE(control$adaptEta))
      .adviAdaptEta(.mu0, .omega0, .theta0, .logPopOmega0, .muRefIdx, .thFix, .omFix,
                    control, .seed)
    else if (length(control$etaCandidates) == 1L) control$etaCandidates
    else 0.1

  .res <- adviLoop_(.mu0, .omega0, .theta0, .logPopOmega0, .muRefIdx, .thFix, .omFix,
                    as.integer(control$iters), as.numeric(.seed), .etaScale,
                    as.numeric(control$tau), as.numeric(control$alpha),
                    as.integer(control$nMc), .it0,
                    .sMu, .sOmega, .sTheta, .sLpo)
  .res$etaScale <- .etaScale
  .res$prep <- .prep
  .res$etaNames <- .prep$etaNames
  .res$thetaNames <- names(.prep$th)
  .res$popOmega <- exp(.res$logPopOmega)
  .res$seed <- .seed
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

  ## point-estimate SEs come from the FOCEi covariance step; covMethod="advi"
  ## (the full-Bayes variational covariance) falls back to "r,s" here.
  .covM <- if (identical(.control$covMethod, "advi")) "r,s" else .control$covMethod
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
  nmObjHandleControlObject(.control, .ret)   # store adviControl for nmObjGetControl.advi
  .ret$foceiModel <- .ui2$focei
  .fit <- nlmixr2CreateOutputFromUi(.ret$ui, data = .ret$origData, control = .fc,
                                    table = .ret$table, env = .ret, est = "advi")
  ## ADVI artifacts + warm-resume state on the fit env
  .e <- .fit$env
  .e$adviElbo <- res$elbo
  .e$adviState <- list(mu = res$mu, omega = res$omega, theta = res$theta,
                       logPopOmega = res$logPopOmega, it0 = res$it0,
                       sMu = res$sMu, sOmega = res$sOmega, sTheta = res$sTheta,
                       sLpo = res$sLpo, seed = res$seed, etaScale = res$etaScale)
  .fit
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
