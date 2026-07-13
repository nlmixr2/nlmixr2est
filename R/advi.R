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
#' @param runLoop function(eta, nAdapt) -> a loop result (list with $elbo,
#'   $parHist), from the initial state; family-agnostic.
#' @noRd
.adviAdaptEta <- function(runLoop, cands, nAdapt) {
  if (length(cands) == 1L) return(cands)
  .best <- -Inf; .bestEta <- cands[1]
  for (.e in cands) {
    .r <- runLoop(.e, nAdapt)
    .el <- .r$elbo
    ## reject a candidate that diverges (non-finite, or the population estimates
    ## in parHist blow up) even if its early ELBO looked good
    .diverged <- any(!is.finite(.el)) || any(!is.finite(.r$parHist)) ||
      max(abs(.r$parHist)) > 1e4
    .score <- if (.diverged) -Inf
      else mean(.el[max(1L, length(.el) - nAdapt %/% 3L):length(.el)])
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
  ## a resumed run keeps its original family; otherwise use the control's
  .fr <- if (!is.null(resume) && !is.null(resume$family))
    identical(resume$family, "fullRank") else identical(control$adviFamily, "fullRank")
  .nL <- neta * (neta + 1L) / 2L

  ## common population/mean init (or resume)
  if (is.null(resume)) {
    .logPopOmega0 <- log(.prep$omega); .mu0 <- matrix(0, N, neta)
    .theta0 <- .prep$theta; .it0 <- 0L
    .sMu <- matrix(0, N, neta); .sTheta <- numeric(ntheta); .sLpo <- numeric(neta)
    if (.fr) {
      ## full-rank: L_i starts diagonal with L_kk = sqrt(popOmega_k) (packed)
      .scale0 <- matrix(0, N, .nL)
      for (.k in seq_len(neta)) .scale0[, .k * (.k + 1L) / 2L] <- exp(0.5 * .logPopOmega0[.k])
      .sScale <- matrix(0, N, .nL)
    } else {
      ## mean-field: log-sd starts at the prior scale
      .scale0 <- matrix(rep(0.5 * .logPopOmega0, each = N), N, neta)
      .sScale <- matrix(0, N, neta)
    }
  } else {
    .mu0 <- resume$mu; .theta0 <- resume$theta; .logPopOmega0 <- resume$logPopOmega
    .it0 <- as.integer(resume$it0); .sMu <- resume$sMu
    .sTheta <- resume$sTheta; .sLpo <- resume$sLpo
    ## full-Bayes stores a generic $scale; point-estimate stores $Lpack/$omega
    .scale0 <- if (!is.null(resume$scale)) resume$scale
      else if (.fr) resume$Lpack else resume$omega
    .sScale <- resume$sScale
  }

  .adviInnerSetup(ui, data, .mu0, control)
  on.exit(.adviInnerFree(), add = TRUE)

  ## the counter-based RNG is keyed by the global iteration index, so resuming
  ## with the original seed continues the exact same stream (prefix property).
  .seed <- if (!is.null(resume) && !is.null(resume$seed)) resume$seed else control$seed
  .muRefIdx <- as.integer(.prep$muRefThetaIdx)
  .thFix <- as.logical(.prep$thetaFix); .omFix <- as.logical(.prep$omegaFix)
  ## per-theta recentering eta (0-based) for mu-referenced intercepts, else -1
  .thetaMuRefEta <- rep(-1L, ntheta)
  for (.k in seq_len(neta)) {
    .p <- .prep$muRefThetaIdx[.k]
    if (!is.na(.p)) .thetaMuRefEta[.p] <- .k - 1L
  }

  ## ---- full-Bayes: variational posterior over the free population vector ----
  if (!isTRUE(control$pointEstimate)) {
    .thetaFreeIdx <- which(!.thFix)               # 1-based ntheta (estimated thetas)
    .omFreeIdx <- which(!.omFix)                   # 1-based eta (estimated variances)
    .npop <- length(.thetaFreeIdx) + length(.omFreeIdx)
    ## phi = c(theta[free], logPopOmega[free]); component -> (theta|omega) maps (0-based)
    .phiThetaIdx <- c(.thetaFreeIdx - 1L, rep(-1L, length(.omFreeIdx)))
    .phiOmIdx <- c(rep(-1L, length(.thetaFreeIdx)), .omFreeIdx - 1L)
    ## phiMuRef[j] = 0-based eta recentered by phi j when it is a mu-ref theta, else -1
    .phiMuRef <- rep(-1L, .npop)
    for (.j in seq_along(.thetaFreeIdx)) {
      .kk <- .thetaMuRefEta[.thetaFreeIdx[.j]]
      if (.kk >= 0) .phiMuRef[.j] <- .kk
    }
    if (is.null(resume)) {
      .mPop0 <- c(.theta0[.thetaFreeIdx], .logPopOmega0[.omFreeIdx])
      .nLpop <- .npop * (.npop + 1L) / 2L
      .Lpop0 <- numeric(.nLpop)
      for (.k in seq_len(.npop)) .Lpop0[.k * (.k + 1L) / 2L] <- 0.1   # init pop posterior sd
      .smPop <- numeric(.npop); .sLpop <- numeric(.nLpop)
    } else {
      .mPop0 <- resume$mPop; .Lpop0 <- resume$Lpop
      .smPop <- resume$smPop; .sLpop <- resume$sLpop
    }
    .runFB <- function(eta, iters, it0 = 0L) {
      adviLoopFB_(.mu0, .scale0, .theta0, .logPopOmega0, .mPop0, .Lpop0,
                  as.integer(.phiThetaIdx), as.integer(.phiOmIdx), as.integer(.phiMuRef),
                  .muRefIdx, as.integer(.fr),
                  as.integer(iters), as.numeric(.seed), eta, as.numeric(control$tau),
                  as.numeric(control$alpha), as.integer(control$nMc), it0,
                  .sMu, .sScale, .smPop, .sLpop)
    }
    .etaScale <- if (!is.null(resume) && !is.null(resume$etaScale)) resume$etaScale
      else if (isTRUE(control$adaptEta))
        .adviAdaptEta(function(e, n) .runFB(e, n), control$etaCandidates,
                      as.integer(min(control$iters, 75L)))
      else if (length(control$etaCandidates) == 1L) control$etaCandidates else 0.1
    .res <- .runFB(.etaScale, control$iters, it0 = .it0)
    .res$family <- control$adviFamily
    .res$pointEstimate <- FALSE
    .res$etaScale <- .etaScale
    .res$prep <- .prep
    .res$etaNames <- .prep$etaNames
    .res$thetaNames <- names(.prep$th)
    .res$popOmega <- exp(.res$logPopOmega)
    .res$seed <- .seed
    .res$phiThetaIdx <- .phiThetaIdx; .res$phiOmIdx <- .phiOmIdx
    ## population variational covariance in phi space (Lpop Lpop^T)
    .Lp <- matrix(0, .npop, .npop)
    for (.i in seq_len(.npop)) for (.j in seq_len(.i)) .Lp[.i, .j] <- .res$Lpop[.i * (.i - 1L) / 2L + .j]
    .res$adviCov <- .Lp %*% t(.Lp)
    class(.res) <- "nlmixr2advi"
    return(.res)
  }

  ## family-appropriate one-loop runner (from the initial state)
  .runLoop <- function(eta, iters, it0 = 0L, scale = .scale0, mu = .mu0, theta = .theta0,
                       lpo = .logPopOmega0, sMu = .sMu, sScale = .sScale,
                       sTheta = .sTheta, sLpo = .sLpo) {
    .fn <- if (.fr) adviLoopFR_ else adviLoop_
    .fn(mu, scale, theta, lpo, .muRefIdx, .thetaMuRefEta, .thFix, .omFix,
        as.integer(iters), as.numeric(.seed), eta, as.numeric(control$tau),
        as.numeric(control$alpha), as.integer(control$nMc), it0,
        sMu, sScale, sTheta, sLpo)
  }
  ## step-size scale: reuse the resumed run's, else adaptively search, else fixed.
  .etaScale <- if (!is.null(resume) && !is.null(resume$etaScale)) resume$etaScale
    else if (isTRUE(control$adaptEta))
      .adviAdaptEta(function(e, n) .runLoop(e, n), control$etaCandidates,
                    as.integer(min(control$iters, 75L)))
    else if (length(control$etaCandidates) == 1L) control$etaCandidates
    else 0.1

  .res <- .runLoop(.etaScale, control$iters, it0 = .it0)
  ## normalize the per-subject scale field name (Lpack for full-rank, omega else)
  .res$scale <- if (.fr) .res$Lpack else .res$omega
  .res$sScale <- .res$sL; if (is.null(.res$sScale)) .res$sScale <- .res$sOmega
  .res$family <- control$adviFamily
  .res$pointEstimate <- TRUE
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
    .nm <- c(res$thetaNames[res$phiThetaIdx[res$phiThetaIdx >= 0] + 1L],
             paste0("omega.", res$etaNames[res$phiOmIdx[res$phiOmIdx >= 0] + 1L]))
    if (nrow(.cov) == length(.nm)) dimnames(.cov) <- list(.nm, .nm)
    .e$adviCov <- .cov
  }
  .e$adviState <- .st
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
