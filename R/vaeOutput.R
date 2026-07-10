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
      enc <- if (prep$covType[j] == "continuous")
        paste0("log(", covNames[j], "/", center, ")")
      else paste0("(", covNames[j], " - ", center, ")")
      terms <- c(terms, paste0(bn, " * ", enc))
      betaVals[[bn]] <- fit$beta[k, j]
    }
    repl <- paste0("(", thName, " + ", paste(terms, collapse = " + "), ")")
    newTxt <- gsub(paste0("\\b", thName, "\\b"), repl, deparse1(.lines[[.idx]]))
    ui2 <- do.call(rxode2::model, list(ui2, str2lang(newTxt)))
  }

  ## 2. set ini() estimates to the VAE solution
  .setIni <- function(u, expr) do.call(rxode2::ini, list(u, str2lang(expr)))
  for (k in seq_along(thetaNames)) {
    ui2 <- .setIni(ui2, paste0(thetaNames[k], " <- ", signif(fit$zPop[k], 12)))
    ui2 <- .setIni(ui2, paste0(fit$prep$etaNames[k], " ~ ", signif(fit$omega[k], 12)))
  }
  for (bn in names(betaVals)) ui2 <- .setIni(ui2, paste0(bn, " <- ", signif(betaVals[[bn]], 12)))
  .errRow <- ui$iniDf[!is.na(ui$iniDf$err) & !is.na(ui$iniDf$ntheta), , drop = FALSE]
  for (en in .errRow$name) {
    .v <- if (!is.null(names(fit$a)) && en %in% names(fit$a)) fit$a[[en]] else fit$a[1]
    ui2 <- .setIni(ui2, paste0(en, " <- ", signif(.v, 12)))
  }
  ui2
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
  .ui2 <- .vaeUpdateModel(.ui, fit)
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
  ## the ORIGINAL (pre-covariate-selection) model for $uiIni/$iniDf0; must be set
  ## AFTER assembly (.nlmixr2FitUpdateParams overwrites $iniDf0 with the global
  ## iniDf data.frame, which cannot represent the structure change)
  .e <- .fit$env
  .e$iniDf0 <- rxode2::rxUiCompress(rxode2::rxUiDecompress(.ui))
  .fit
}
