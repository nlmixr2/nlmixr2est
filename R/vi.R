# vi.R -- orchestration for est="emvi" (variational EM) and est="fbvi" (full
# Bayes VI), in the style of Kucukelbir et al. 2017 but not the published
# algorithm (see emviControl()).  Sets up the FOCEi inner problem (reused for the
# per-subject log-joint and eta-gradient) plus, when non-mu structural thetas are
# present, the impmap theta-sensitivity model (reused for the outer population
# gradient), then drives the optimization loop in C++.
#
# The internal names keep the historical `advi` spelling -- the C++ entry point
# is still adviOptimize_ and the FOCEi inner marker is still est="advi" -- so
# only the user-facing surface moved.  The inner marker in particular is load
# bearing: .foceiOptEnvLik selects the theta-sensitivity model build on it.

#' A foceiControl carrying the chosen inner likelihood + solving options.
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
#' @param control emviControl
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
  ## declared population omega structure: installs the off-diagonal mask so the
  ## ELBO/gradient entry points see it without going through adviOptimize_
  .ob <- .omegaBlockFromIniDf(.ui$iniDf, .foceiEtaThetaMap(.ui)$etaNames)
  .env$adviOmegaMat <- .ob$mat
  .env$adviOmegaFixMat <- .ob$fixMat
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

#' Run the variational optimization: prep, inner setup, initialize the variational +
#' population state, and drive the whole optimization (the adaptive step-size
#' search + the main loop) in one C++ call (adviOptimize_).
#' @param ui bounded-transformed rxode2 ui
#' @param data estimation data
#' @param control emviControl
#' @param resume optional list from a previous fit's `$viState` for warm resume
#' @return the raw result list (variational params, estimates, elbo, parHist)
#' @noRd
.adviOptimize <- function(ui, data, control, resume = NULL) {
  .prep <- .adviDataPrep(ui, data)
  N <- .prep$N; neta <- .prep$neta
  ## a resumed run keeps its original family; otherwise use the control's
  .fr <- if (!is.null(resume) && !is.null(resume$family))
    identical(resume$family, "fullRank") else identical(control$viFamily, "fullRank")

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
  ## pointEstimate picks the ALGORITHM, and every caller resolves it from `est`
  ## before getting here.  An unresolved NULL would slip through isTRUE() as
  ## FALSE and quietly run full Bayes under est="emvi", so refuse instead.
  if (!isTRUE(control$pointEstimate) && !isFALSE(control$pointEstimate)) {
    stop("emviControl(pointEstimate=) was never resolved from `est`", call. = FALSE)
  }
  ## correlated etas: the point-estimate families estimate the full omega block;
  ## the full-Bayes path (pointEstimate=FALSE) parameterizes phi with per-eta
  ## log-variances only, so it cannot carry an off-diagonal yet
  if (.omegaHasOffDiag(.prep$omegaMat) && !isTRUE(control$pointEstimate)) {
    stop("est=\"fbvi\" does not support correlated etas", call. = FALSE)
  }
  .res <- adviOptimize_(list(
    pointEstimate = isTRUE(control$pointEstimate), fr = as.integer(.fr),
    N = as.integer(N),
    theta = as.numeric(.prep$theta), omega = as.numeric(.prep$omega),
    omegaMat = .prep$omegaMat, omegaFixMat = .prep$omegaFixMat,
    perNoCor = as.numeric(control$perNoCor),
    tol = as.numeric(control$tol), evalElbo = as.integer(control$evalElbo),
    jacType = as.integer(.prep$jacType), jacRange = as.numeric(.prep$jacRange),
    klWarmup = as.integer(control$klWarmup), temperInit = as.numeric(control$temperInit),
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
  .res$family <- control$viFamily
  .res$prep <- .prep
  .res$etaNames <- .prep$etaNames
  .res$thetaNames <- names(.prep$th)
  .res$model <- .setup$model
  class(.res) <- "nlmixr2vi"
  .res
}

#' Assemble the standard nlmixr2FitData from a variational result: seed the ui
#' iniDf with the estimates (population thetas + between-subject omega), supply
#' the variational posterior means as the FOCEi inner EBE start (etaMat), and run
#' the eval-only FOCEi finalize (maxOuterIterations=0) which reuses inner.cpp for
#' the objective, EBEs, residual tables, and the covariance step.  No outer
#' optimizer is run; the variational estimates are final.  Mirrors .vaeToFit.
#' @noRd
.adviToFit <- function(env, res) {
  .ui <- env$ui
  .control <- env$emviControl
  ## which of the two methods produced this fit; env$est is set by .viEst, but
  ## fall back to pointEstimate so a directly-called .adviFitModel still labels
  ## the fit with a real method name rather than NULL
  .est <- env$est
  if (!is.character(.est) || length(.est) != 1L || !(.est %in% c("emvi", "fbvi"))) {
    .est <- if (isTRUE(res$pointEstimate)) "emvi" else "fbvi"
  }
  .prep <- res$prep
  .rxControl <- .control$rxControl

  .ret <- new.env(parent = emptyenv())
  .ret$table <- env$table
  .foceiPreProcessData(env$data, .ret, .ui, .rxControl)

  ## seed the ui iniDf with the variational estimates so the eval reports them
  .uiD <- rxode2::rxUiDecompress(.ui)
  .idf <- .uiD$iniDf
  .thRow <- !is.na(.idf$ntheta)
  .idf$est[.thRow] <- res$theta[.idf$ntheta[.thRow]]
  .popOm <- stats::setNames(res$popOmega, .prep$etaNames)
  .etaRow <- !is.na(.idf$neta1) & .idf$neta1 == .idf$neta2
  .idf$est[.etaRow] <- .popOm[.idf$name[.etaRow]]
  ## estimated off-diagonals (full omega block), keyed by the iniDf neta indices
  .omM <- if (is.null(res$popOmegaMat)) diag(res$popOmega, .prep$neta) else res$popOmegaMat
  dimnames(.omM) <- list(.prep$etaNames, .prep$etaNames)
  .diagRow <- .idf[.etaRow, , drop = FALSE]
  .etaIdx <- stats::setNames(match(.diagRow$name, .prep$etaNames),
                             as.character(.diagRow$neta1))
  .offRow <- which(!is.na(.idf$neta1) & .idf$neta1 != .idf$neta2)
  for (.r in .offRow) {
    .i <- .etaIdx[as.character(.idf$neta1[.r])]
    .j <- .etaIdx[as.character(.idf$neta2[.r])]
    if (!is.na(.i) && !is.na(.j)) .idf$est[.r] <- .omM[.i, .j]
  }
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
  .ret$omega <- .omM
  .ret$ui <- .ui2
  .ret$fullTheta <- stats::setNames(res$theta, names(.prep$th))

  ## covMethod="vi": for full-Bayes the SEs come from the population variational
  ## covariance (installed below, so skip the FOCEi cov step); for point-estimate
  ## there is no population variational block, so fall back to the FOCEi "r,s".
  .covM <- if (identical(.control$covMethod, "vi"))
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
  .ret$method <- .est
  .ret$extra <- ""
  .ret$est <- .est
  .ret$ofvType <- .est
  .ret$adjObf <- .control$adjObf
  ## the optimization walk (standard parHistData -> $parHist accessor)
  if (!is.null(res$parHistData)) .ret$parHistData <- res$parHistData
  nmObjHandleControlObject(.control, .ret)   # store emviControl for nmObjGetControl.advi
  ## reuse the models compiled for the variational loop (inner/EBE/pred + thetaSens):
  ## with $model present the eval-only finalize skips its own symengine rebuild
  ## (the finalize reads only the foce-prefix columns of the inner model, so the
  ## interaction-model column layout is compatible)
  if (!is.null(res$model)) {
    .ret$model <- res$model
  } else {
    .ret$foceiModel <- .ui2$focei
  }
  .fit <- nlmixr2CreateOutputFromUi(.ret$ui, data = .ret$origData, control = .fc,
                                    table = .ret$table, env = .ret, est = .est)
  ## variational artifacts + warm-resume state on the fit env
  .e <- .fit$env
  .e$viElbo <- res$elbo
  ## an early ELBO-convergence stop is a real difference from the requested
  ## `iters`; say so rather than leaving a short trace to be noticed
  if (isTRUE(res$tolStopped)) {
    warning(sprintf("ELBO converged at iteration %d of %d (emviControl(tol=))",
                    length(res$elbo), as.integer(.control$iters)),
            call. = FALSE)
  }
  ## The adaptEta search picking the largest (or smallest) candidate is the one
  ## case where the grid itself is plausibly the binding constraint -- the model
  ## may want a step outside it and cannot say so.  Surface that instead of
  ## letting it look like a converged choice; $etaScores shows the full search.
  .cand <- as.numeric(.control$etaCandidates)
  if (length(res$etaScores) > 1L && length(.cand) > 1L &&
        any(is.finite(res$etaScores))) {
    ## EXACT comparison, not all.equal: etaScale is assigned straight from an
    ## element of etaCandidates in C++, so it is bit-identical, and all.equal's
    ## relative tolerance would call two genuinely distinct neighbouring
    ## candidates equal -- c(0.1, 0.1 + 1e-9) selecting the bottom would be
    ## reported as the top.
    .sel <- as.numeric(res$etaScale)
    if (isTRUE(.sel == max(.cand))) {
      warning("step-size search hit the top of etaCandidates; consider widening",
              call. = FALSE)
    } else if (isTRUE(.sel == min(.cand))) {
      warning("step-size search hit the bottom of etaCandidates; consider widening",
              call. = FALSE)
    }
  }
  .st <- list(mu = res$mu, theta = res$theta, logPopOmega = res$logPopOmega,
              popOmegaMat = res$popOmegaMat, nbCorrel = res$nbCorrel,
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
    .cov <- res$viCov
    ## adviOptimize_ always returns viCov on the full-Bayes branch; a missing one
    ## means the C++ result contract and this reader have drifted apart, which
    ## otherwise surfaces as `nrow(NULL)` -> "argument is of length zero"
    if (!is.matrix(.cov)) {
      stop("full Bayes returned no variational covariance ('viCov')", call. = FALSE)
    }
    .thNm <- res$prep$thetaRealNames[res$phiThetaIdx[res$phiThetaIdx >= 0] + 1L]
    .nm <- c(.thNm, paste0("omega.", res$etaNames[res$phiOmIdx[res$phiOmIdx >= 0] + 1L]))
    if (nrow(.cov) == length(.nm)) dimnames(.cov) <- list(.nm, .nm)
    .e$viCov <- .cov
    ## covMethod="vi": install the population variational covariance as the
    ## fit's SE source (the theta block maps directly to parFixedDf's
    ## population/residual parameters).  For any other covMethod the FOCEi
    ## covariance step (analytic/r,s/...) already ran on the full inner model;
    ## only fall back to the variational covariance if that chain came up empty.
    .cmDone <- tryCatch(as.character(.e$covMethod), error = function(e) "")
    if (identical(.control$covMethod, "vi")) {
      .adviInstallVarCov(.fit, res)
    } else if (is.null(.e$cov) || !is.matrix(.e$cov) ||
                 length(.cmDone) != 1L || !nzchar(.cmDone) ||
                 identical(.cmDone, "failed")) {
      message("covMethod=\"", .control$covMethod,
              "\" covariance was not available; using the variational covariance")
      .adviInstallVarCov(.fit, res)
    }
  }
  .e$viState <- .st
  .fit
}

#' Install a named theta-block covariance as the fit's covariance and update
#' the parFixedDf SEs / RSEs / CIs accordingly.
#' @noRd
.adviInstallThetaCov <- function(fit, cov) {
  .env <- if (rxode2::rxIs(fit, "nlmixr2FitData")) fit$env else fit
  if (!is.environment(.env) || is.null(cov)) return(invisible())
  .env$cov <- cov
  if (!exists("parFixedDf", envir = .env, inherits = FALSE)) return(invisible())
  .pf <- .env$parFixedDf
  .se <- sqrt(diag(cov))
  .ci <- tryCatch(as.numeric(rxode2::rxGetControl(.env$ui, "ci", 0.95)), error = function(e) 0.95)
  .qn <- stats::qnorm(1 - (1 - .ci) / 2)
  for (.n in rownames(.pf)) {
    if (.n %in% names(.se) && "SE" %in% names(.pf)) {
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
  invisible()
}

#' Install the population variational covariance (Lpop Lpop^T) as the fit's
#' covariance + parFixedDf SEs -- the natural full-Bayes uncertainty.  The theta
#' block of the phi-space covariance maps directly to the population / residual
#' parameters (by name); the log-variance block is retained on $env$viCov.
#' @noRd
.adviInstallVarCov <- function(fit, res) {
  .thComp <- which(res$phiThetaIdx >= 0)
  if (length(.thComp) == 0L) return(invisible())
  .thNames <- res$prep$thetaRealNames[res$phiThetaIdx[.thComp] + 1L]
  .thetaCov <- res$viCov[.thComp, .thComp, drop = FALSE]
  dimnames(.thetaCov) <- list(.thNames, .thNames)
  .adviInstallThetaCov(fit, .thetaCov)
  .env <- if (rxode2::rxIs(fit, "nlmixr2FitData")) fit$env else fit
  if (is.environment(.env)) .env$covMethod <- "vi"
  invisible()
}

#' Fit an emvi/fbvi model: set up the inner/outer problems and run the C++ loop.
#' @param env estimation environment (holds ui, data, emviControl)
#' @noRd
.adviFitModel <- function(env) {
  .ui <- env$ui
  .control <- env$emviControl
  ## warm resume: accept a prior emvi/fbvi fit or its viState
  .resume <- .control$resume
  if (!is.null(.resume)) {
    if (rxode2::rxIs(.resume, "nlmixr2FitData")) .resume <- .resume$env$viState
    else if (is.environment(.resume) && exists("viState", .resume)) .resume <- .resume$viState
    if (!is.list(.resume) || is.null(.resume$it0))
      stop("'resume' must be a prior emvi/fbvi fit or its $env$viState", call. = FALSE)
    ## the two methods carry DIFFERENT state: emvi saves sTheta/sLpo, fbvi saves
    ## mPop/Lpop/smPop/sLpop.  Resuming across them reaches adviOptimize_ with the
    ## wrong half missing and dies on an Rcpp NULL conversion, so say what is
    ## actually wrong instead.
    if (!is.null(.resume$pointEstimate) &&
          !identical(isTRUE(.resume$pointEstimate), isTRUE(.control$pointEstimate))) {
      stop("cannot resume an ", if (isTRUE(.resume$pointEstimate)) "emvi" else "fbvi",
           " fit under est=\"", if (isTRUE(.control$pointEstimate)) "emvi" else "fbvi",
           "\"", call. = FALSE)
    }
  }
  .res <- .adviOptimize(.ui, env$data, .control, resume = .resume)
  if (isTRUE(.control$returnVi)) return(.res)
  .adviToFit(env, .res)
}
