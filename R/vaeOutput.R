# vaeOutput.R -- build the fitted nlmixr2 output from a trained VAE. The model is
# updated with the selected covariate effects (exact centered expressions, e.g.
# ka <- exp(lka + beta_lka_WT*log(WT/79.6) + eta.ka)) and the ini() estimates are
# set to the VAE solution. The standard nlmixr2FitData (objective, parFixed/SEs,
# EBEs, residuals, tables) is then assembled with nlmixr2CreateOutputFromUi
# driving the FOCEi inner problem at the VAE estimates (no outer optimizer),
# which reuses inner.cpp's engine wholesale -- including the focei covariance
# step ("analytic", "r,s", "r", "s") -- see .vaeToFit.

#' Update a ui with the VAE's selected covariate effects and fitted estimates.
#'
#' Continuous covariates enter as `beta*log(COV/center)`, categorical as
#' `beta*(COV - center)`, inserted into each mu-referenced parameter's model
#' line. Returns the updated ui.
#' @noRd
.vaeUpdateModel <- function(ui, fit) {
  prep <- fit$prep
  .map <- .foceiEtaThetaMap(ui)
  thetaNames <- .map$thetaForEta        # mu-referenced theta per eta (e.g. lka)
  covNames <- fit$covNames
  ui2 <- ui
  betaVals <- list()

  ## 1. inject covariate terms into each parameter's model line
  for (k in seq_along(thetaNames)) {
    ## a free/fixed eta (thetaForEta == NA: literalFix-ed or non-mu-referenced)
    ## has no structural theta to attach a covariate to and is excluded from
    ## covariate selection -- skip it (its structure is already in the model)
    if (is.na(thetaNames[k])) next
    sel <- if (is.null(fit$selected)) integer(0) else which(fit$selected[k, ])
    if (length(sel) == 0L) next
    thName <- thetaNames[k]
    .lines <- ui2$lstExpr
    .idx <- which(vapply(.lines, function(e) thName %in% all.vars(e), logical(1)))
    if (length(.idx) == 0L) next
    .idx <- .idx[1]
    terms <- character(0)
    for (j in sel) {
      bn <- paste0("beta_", thName, "_", covNames[j])
      center <- signif(prep$covPop[j], 12)
      ## an uncentered covariate (a 0/1 indicator) enters bare -- no "- 0" term
      enc <- if (prep$covType[j] == "continuous") {
        paste0("log(", covNames[j], "/", center, ")")
      } else if (center == 0) {
        covNames[j]
      } else {
        paste0("(", covNames[j], " - ", center, ")")
      }
      terms <- c(terms, paste0(bn, " * ", enc))
      betaVals[[bn]] <- fit$beta[k, j]
    }
    ## Inject the covariate terms FLAT (no wrapping parentheses): the mu-ref line
    ## is `p <- exp(theta + eta)`, so replacing `theta` with `theta + beta*cov`
    ## keeps the additive `exp(theta + beta*cov + eta)` form rxode2 recognizes as
    ## a mu-referenced exp() parameter.  Wrapping in parens -- `exp((theta +
    ## beta*cov) + eta)` -- hides the exp() back-transform from muRefCurEval, so
    ## the theta prints on the raw log scale instead of back-transformed.
    repl <- paste0(thName, " + ", paste(terms, collapse = " + "))
    newTxt <- gsub(paste0("\\b", thName, "\\b"), repl, deparse1(.lines[[.idx]]))
    ui2 <- do.call(rxode2::model, list(ui2, str2lang(newTxt)))
  }

  ## 2. set ini() estimates to the VAE solution
  .setIni <- function(u, expr) do.call(rxode2::ini, list(u, str2lang(expr)))
  for (k in seq_along(thetaNames)) {
    ## a free/fixed eta has no structural theta (thetaForEta == NA) -- its
    ## population location is already a literal in the model, so only set omega
    if (!is.na(thetaNames[k])) {
      ui2 <- .setIni(ui2, paste0(thetaNames[k], " <- ", signif(fit$zPop[k], 12)))
    }
    ui2 <- .setIni(ui2, paste0(fit$prep$etaNames[k], " ~ ", signif(fit$omega[k], 12)))
  }
  for (bn in names(betaVals)) ui2 <- .setIni(ui2, paste0(bn, " <- ", signif(betaVals[[bn]], 12)))
  .errRow <- ui$iniDf[!is.na(ui$iniDf$err) & !is.na(ui$iniDf$ntheta), , drop = FALSE]
  for (en in .errRow$name) {
    .v <- if (!is.null(names(fit$a)) && en %in% names(fit$a)) fit$a[[en]] else fit$a[1]
    ui2 <- .setIni(ui2, paste0(en, " <- ", signif(.v, 12)))
  }
  ## 3. non-mu thetas estimated by the bobyqa regression (nonMuTheta="regress"):
  ## these have no eta, so write each regressed value straight into its ini() est.
  if (!is.null(fit$regressTheta) && length(fit$regressTheta) > 0L &&
      !is.null(names(fit$regressTheta))) {
    for (rn in names(fit$regressTheta)) {
      .rv <- fit$regressTheta[[rn]]
      if (is.finite(.rv)) ui2 <- .setIni(ui2, paste0(rn, " <- ", signif(.rv, 12)))
    }
  }
  ## The incremental model()/ini() edits above leave the ui's cached `covariates`
  ## stale: an injected covariate-coefficient theta (beta_<par>_<cov>) is added to
  ## the iniDf as a theta but ALSO stays listed as a covariate.  The augmented
  ## covariance solve then declares that beta_ both as its THETA[k] and as a
  ## phantom data covariate, and fails ("required for solving: beta_...").
  ## Re-parsing the accumulated model function yields a consistent theta/covariate
  ## classification (verified: only the true data covariates remain).
  rxode2::assertRxUi(ui2$fun)
}

#' Update a PINNED VAE fit's model with the estimates.
#'
#' Unlike `.vaeUpdateModel` (which injects fresh `beta_<par>_<cov>` terms into a
#' covariate-free base model), the pinned path keeps the user's ORIGINAL model
#' -- their covariate terms, coefficient names and centers stay exactly as
#' written -- and only writes ini() estimates.  A declared covariate the search
#' selected gets its estimated slope; one it dropped is set to `0` (the term
#' stays in the model).  Pinned covariates are searched at their MODEL value (the
#' model's own centering is retained, e.g. from mu2/mu3 `nlmixrMuDerCov#`, with no
#' extra VAE mean-centering), so `zPop` is already the model intercept and the
#' coefficient transfers directly with no correction.  Out-of-pool declared
#' covariates were estimated in place by the regress M-step (`regressTheta`).
#' @noRd
.vaeUpdateModelPinned <- function(ui, fit) {
  prep <- fit$prep
  pairs <- prep$pinPairs
  thetaNames <- .foceiEtaThetaMap(ui)$thetaForEta   # mu-referenced theta per eta
  covNames <- fit$covNames
  ui2 <- ui
  .setIni <- function(u, expr) do.call(rxode2::ini, list(u, str2lang(expr)))

  ## 1. in-pool declared coefficients: selected -> estimated slope, dropped -> 0.
  ## The pinned covariates are searched at their MODEL value (no VAE re-centering,
  ## the model's own centering is retained), so zPop is already the model
  ## intercept -- the coefficient transfers directly with no correction.
  .inRows <- if (is.null(pairs)) NULL else pairs[pairs$inPool, , drop = FALSE]
  for (.r in seq_len(NROW(.inRows))) {
    .k <- .inRows$k[.r]
    .j <- match(.inRows$covName[.r], covNames)
    .sel <- !is.null(fit$selected) && !is.na(.j) && isTRUE(fit$selected[.k, .j])
    .betaVal <- if (.sel) fit$beta[.k, .j] else 0
    ui2 <- .setIni(ui2, paste0(.inRows$coefName[.r], " <- ", signif(.betaVal, 12)))
  }

  ## 2. structural population thetas + omega per eta
  for (k in seq_along(thetaNames)) {
    if (!is.na(thetaNames[k])) {
      ui2 <- .setIni(ui2, paste0(thetaNames[k], " <- ", signif(fit$zPop[k], 12)))
    }
    ui2 <- .setIni(ui2, paste0(prep$etaNames[k], " ~ ", signif(fit$omega[k], 12)))
  }

  ## 3. residual error params
  .errRow <- ui$iniDf[!is.na(ui$iniDf$err) & !is.na(ui$iniDf$ntheta), , drop = FALSE]
  for (en in .errRow$name) {
    .v <- if (!is.null(names(fit$a)) && en %in% names(fit$a)) fit$a[[en]] else fit$a[1]
    ui2 <- .setIni(ui2, paste0(en, " <- ", signif(.v, 12)))
  }

  ## 4. out-of-pool declared covariates + non-mu thetas estimated by the regress
  ## M-step (written straight into their ini() est)
  if (!is.null(fit$regressTheta) && length(fit$regressTheta) > 0L &&
        !is.null(names(fit$regressTheta))) {
    for (rn in names(fit$regressTheta)) {
      .rv <- fit$regressTheta[[rn]]
      if (is.finite(.rv)) ui2 <- .setIni(ui2, paste0(rn, " <- ", signif(.rv, 12)))
    }
  }
  rxode2::assertRxUi(ui2$fun)
}

#' Translate the vaeControl into the foceiControl that drives the output step:
#' no outer/inner optimization (the VAE estimates and encoder etas are final),
#' the VAE's chosen inner likelihood (focei -> interaction=1; foce/focep ->
#' interaction=0, focep = FOCE+ with R at the live conditional eta), and the
#' covMethod passed through so the focei covariance step ("analytic", "r,s",
#' "r", "s", "") runs directly on the frozen problem.
#' @noRd
.vaeControlToFoceiControl <- function(env, assign = TRUE) {
  .control <- env$vaeControl
  .lik <- .control$likelihood
  .interaction <- if (.lik %in% c("foce", "focep")) 0L else 1L
  .foce <- if (identical(.lik, "focep")) "foce+" else "nonmem"
  .fc <- foceiControl(rxControl = .control$rxControl,
                      maxOuterIterations = 0L, maxInnerIterations = 0L,
                      covMethod = .control$covMethod,
                      etaMat = env$etaMat,
                      interaction = .interaction, foce = .foce,
                      sumProd = .control$sumProd,
                      optExpression = .control$optExpression,
                      literalFix = .control$literalFix,
                      literalFixRes = .control$literalFixRes,
                      addProp = .control$addProp,
                      calcTables = .control$calcTables,
                      compress = .control$compress,
                      ci = .control$ci,
                      sigdigTable = .control$sigdigTable,
                      stickyRecalcN = .control$stickyRecalcN,
                      maxOdeRecalc = .control$maxOdeRecalc,
                      odeRecalcFactor = .control$odeRecalcFactor,
                      indTolRelax = .control$indTolRelax,
                      eventSens = .control$eventSens,
                      fast = FALSE, # no outer optimizer -- skip the outer gradient model
                      print = 0L)
  if (assign) env$control <- .fc
  .fc
}

#' Assemble the standard nlmixr2FitData from a trained VAE with
#' nlmixr2CreateOutputFromUi (the nlme/nlm/nlmer output pattern), driving the
#' FOCEi INNER problem at the VAE's fixed population estimates
#' (maxOuterIterations=0 -- no outer optimizer is run) with the encoder etas
#' supplied as etaMat. This reuses inner.cpp's parallel (OpenMP) inner
#' likelihood wholesale: multiple endpoints, multiple error structures,
#' log-likelihood, M2/M3/M4 censoring, and MIXTURE hard-assignment (nSub*nMix
#' per-component solves via setIndMixest -> mixNum/mixList) -- none of which is
#' reimplemented here. The covariance is computed by the focei covariance step
#' itself (covMethod passed through .vaeControlToFoceiControl) and returned on
#' the fit. The model is first updated with the selected covariate effects; the
#' ORIGINAL (pre-covariate) ui is stashed in $iniDf0 for the iniUi/iniDf0
#' accessors.
#' @noRd
.vaeToFit <- function(env, fit) {
  .ui <- env$ui
  .control <- env$vaeControl
  ## mu2/mu3 restore info staged by the preprocess hook: the focei covariance
  ## recompute below re-runs preprocessing and clears it, so snapshot it here and
  ## reinstate it just before the mu2 finalize restores the original model.  The
  ## on.exit guard keeps the global consistent even if assembly errors out.
  .savedMuRef <- .muRefTrans$cur
  on.exit(.muRefTrans$cur <- .savedMuRef, add = TRUE)
  .ui2 <- if (isTRUE(fit$prep$pinActive)) .vaeUpdateModelPinned(.ui, fit) else .vaeUpdateModel(.ui, fit)
  ## Collapse any etas injected for non-mu-referenced thetas (nonMuTheta="eta"/
  ## "fix"): .vaeUpdateModel has already written the population estimate (zPop =
  ## theta+mean(eta)) into the theta, so drop the temporary eta from the reported
  ## model and its column from the EBE matrix -- the parameter is reported as a
  ## plain fixed effect.
  ## per-fit record (set by the preprocess hook, copied onto env by runPreProcess);
  ## fall back to the global only for direct callers that bypass the hook wrapper
  .injEtas <- if (exists("vaeNonMuEtas", envir = env, inherits = FALSE)) {
    env$vaeNonMuEtas
  } else {
    nlmixr2global$nlmixr2EstEnv$vaeNonMuEtas
  }
  if (length(.injEtas) > 0L) {
    .injEtas <- .injEtas[.injEtas %in% fit$prep$etaNames]
    if (length(.injEtas) > 0L) {
      .ui2 <- rmEta(.ui2, .injEtas)
      .keep <- !(fit$prep$etaNames %in% .injEtas)
      fit$mu <- fit$mu[, .keep, drop = FALSE]
      fit$zPopMat <- fit$zPopMat[, .keep, drop = FALSE]
      fit$prep$etaNames <- fit$prep$etaNames[.keep]
    }
  }
  .ret <- new.env(parent = emptyenv())
  .ret$table <- env$table
  ## encoder etas as the FOCEi inner starting point [nsub, neta] in eta order
  .etaMat <- fit$mu - fit$zPopMat
  colnames(.etaMat) <- fit$prep$etaNames
  .ret$etaMat <- .etaMat
  ## presetting $method/$extra keeps the C++ finalize from writing the focei
  ## "FOCE"/"i (outer: ...)" labels (and from clobbering $parHistData below)
  .ret$method <- "vae"
  .ret$extra <- ""
  .ret$est <- "vae"
  .ret$adjObf <- .control$adjObf
  ## the VAE training artifacts + the ORIGINAL model for $uiIni/$iniDf0
  .ret$vae <- list(elboTrace = fit$elboTrace, beta = fit$beta, selected = fit$selected,
                   covNames = fit$covNames, zPop = fit$zPop, omega = fit$omega, a = fit$a,
                   seed = .control$seed)
  ## the VAE optimization walk (standard parHistData -> $parHist accessor)
  if (!is.null(fit$parHist)) .ret$parHistData <- fit$parHist
  nmObjHandleControlObject(.control, .ret) # stores $vaeControl for nmObjGetControl.vae
  .vaeControlToFoceiControl(.ret)
  .fit <- nlmixr2CreateOutputFromUi(.ui2, data = env$data, control = .ret$control,
                                    table = env$table, env = .ret, est = "vae")
  ## mu2/mu3/mu4 covariate rewriting: restore the original algebraic covariate
  ## expression (nlmixrMuDerCov# -> e.g. wt.cl*(WT/70)) in the reported model and
  ## drop the derived data columns.  VAE assembles its output outside the focei
  ## path that normally runs this, so reinstate the restore info and invoke the
  ## mu2 finalize hook directly.
  .muRefTrans$cur <- .savedMuRef
  .fit <- .uiFinalizeMu2hook(.fit)
  ## the ORIGINAL (pre-covariate-selection) model for $uiIni/$iniDf0; must be set
  ## AFTER assembly (.nlmixr2FitUpdateParams overwrites $iniDf0 with the global
  ## iniDf data.frame, which cannot represent the structure change).  Use the pure
  ## input ui (pre-mu2-rewrite) when available so iniDf0 shows the user's model.
  .e <- .fit$env
  .origUi <- if (!is.null(env$nlmixrPureInputUi)) env$nlmixrPureInputUi else .ui
  .e$iniDf0 <- rxode2::rxUiCompress(rxode2::rxUiDecompress(.origUi))
  .fit
}
