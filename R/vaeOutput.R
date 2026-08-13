# vaeOutput.R -- build the fitted nlmixr2 output from a trained VAE. The model is
# updated with the selected covariate effects (exact centered expressions, e.g.
# ka <- exp(lka + beta.lka.WT.log*log(WT/79.6) + eta.ka)) and the ini() estimates are
# set to the VAE solution. The standard nlmixr2FitData (objective, parFixed/SEs,
# EBEs, residuals, tables) is then assembled with nlmixr2CreateOutputFromUi
# driving the FOCEi inner problem at the VAE estimates (no outer optimizer),
# which reuses inner.cpp's engine wholesale -- including the focei covariance
# step ("analytic", "r,s", "r", "s") -- see .vaeToFit.

#' Would a corrected structural theta still respect its ini() bounds?
#' @param ui rxode2 ui carrying the ini() bounds
#' @param thName structural theta name
#' @param value the corrected estimate
#' @return TRUE when the value is finite and within the declared bounds
#' @noRd
.vaeIcAdjInBounds <- function(ui, thName, value) {
  if (!is.finite(value)) return(FALSE)
  .idf <- tryCatch(ui$iniDf, error = function(e) NULL)
  if (is.null(.idf)) return(TRUE)
  .i <- match(thName, .idf$name)
  if (is.na(.i)) return(TRUE)
  .lo <- if (is.null(.idf$lower)) -Inf else .idf$lower[.i]
  .hi <- if (is.null(.idf$upper)) Inf else .idf$upper[.i]
  if (is.na(.lo)) .lo <- -Inf
  if (is.na(.hi)) .hi <- Inf
  value > .lo && value < .hi
}

#' Make a generated coefficient name unique against those already emitted
#' @noRd
.vaeUniqueName <- function(nm, used) {
  if (!(nm %in% used)) return(nm)
  .i <- 2L
  while (paste0(nm, .i) %in% used) .i <- .i + 1L
  paste0(nm, .i)
}

#' Inject covariate terms beside the mu-referenced structural theta
#'
#' Structural rather than textual.  A `gsub` over the deparsed line replaces
#' EVERY occurrence of the theta, so a model that mentions it again in the same
#' line -- `ka <- exp(lka + eta.ka) + 0 * lka` -- has the covariate terms
#' injected twice, the second time in a position operator precedence turns into
#' an extra additive effect.  The theta is replaced only where it sits in the
#' same additive group as its eta, which is the mu-referenced occurrence.
#' @param line model line (language object)
#' @param thName structural theta name
#' @param etaName that parameter's eta name
#' @param termsTxt covariate terms to add, already collapsed with `+`
#' @return the rewritten line, or `NULL` when no mu-referenced occurrence is found
#' @noRd
.vaeInjectCov <- function(line, thName, etaName, termsTxt) {
  .repl <- str2lang(paste0(thName, " + ", termsTxt))
  .done <- FALSE
  .isAdd <- function(e) is.call(e) && is.name(e[[1L]]) &&
    as.character(e[[1L]]) %in% c("+", "-")
  ## names reachable through a flattened +/- chain
  .addVars <- function(e) {
    if (is.name(e)) return(as.character(e))
    if (.isAdd(e)) return(unlist(lapply(as.list(e)[-1L], .addVars), use.names = FALSE))
    character(0)
  }
  .sub <- function(x) {
    if (.done) return(x)
    if (is.name(x) && identical(as.character(x), thName)) {
      .done <<- TRUE
      return(.repl)
    }
    if (.isAdd(x)) for (.i in seq_along(x)[-1L]) x[[.i]] <- .sub(x[[.i]])
    x
  }
  .rec <- function(e) {
    if (.done || !is.call(e)) return(e)
    if (.isAdd(e)) {
      .v <- .addVars(e)
      if (thName %in% .v && etaName %in% .v) return(.sub(e))
    }
    for (.i in seq_along(e)[-1L]) e[[.i]] <- .rec(e[[.i]])
    e
  }
  .out <- .rec(line)
  if (!.done) NULL else .out
}

#' Which shape to write for a selected covariate column
#'
#' The search picks a FAMILY (column); shapes within a family span the same
#' model, so the written parameterization is the first shape this
#' (parameter, covariate) pair allows in that family, falling back to the
#' column's own representative shape.
#'
#' A hockey arm is exempt, as `cat` is: the family holds exactly one
#' parameterization, and the user-facing name `"hockey"` names the RELATIONSHIP,
#' not either arm -- swapping it in would emit a shape with no expression.
#' @noRd
.vaeSelectedShape <- function(prep, parAliases, j) {
  .own <- prep$covShape[j]
  if (identical(.own, "cat") || .own %in% .vaeHockeyArms ||
        is.null(prep$shapeRules)) return(.own)
  .ok <- .vaeShapesFor(prep$shapeRules, parAliases, prep$covRaw[j])
  ## never write a shape that is not expressible at this center
  .ok <- .ok[.vaeShapeUsable(.ok, prep$covPop[j])]
  if (length(.ok) == 0L) return(.own)
  .m <- .ok[.vaeShapeFamily(.ok) == prep$covFamily[j]]
  ## the covAllow shape mask makes an empty match unreachable in a real fit;
  ## the column's own shape is a valid parameterization of the same family
  if (length(.m) == 0L) .own else .m[1L]
}

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
  ## a shape written in a non-centered parameterization moves part of the effect
  ## into the intercept, so the structural theta is corrected by this much
  icAdj <- rep(0, length(thetaNames))

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
    ## names this parameter answers to in a shapes= rule
    .aliases <- c(prep$etaNames[k], thName, sub("^eta\\.", "", prep$etaNames[k]))
    for (j in sel) {
      ## The selected column fixes the FAMILY; which parameterization of that
      ## family is written is the user's choice, so take the first shape this
      ## (parameter, covariate) pair allows within the family.
      .shp <- .vaeSelectedShape(prep, .aliases, j)
      .ctr <- prep$covPop[j]
      .raw <- prep$covRaw[j]
      .r <- .vaeShapeBeta(.shp, .ctr, fit$beta[k, j])
      ## An uncentered parameterization moves part of the effect into the
      ## intercept, which the C++ M-step has already clamped to the theta's
      ## ini() bounds -- so the corrected value can land outside them.  Clamping
      ## it here would change the prediction, so fall back to the family's
      ## CENTERED shape instead: same fit, no intercept correction needed.
      ## (a hockey arm never gets here: both arms vanish at the knot, so their
      ## interceptAdj is 0 and there is nothing to push out of bounds)
      if (.r$interceptAdj != 0 &&
            !.vaeIcAdjInBounds(ui, thName, fit$zPop[k] + icAdj[k] + .r$interceptAdj)) {
        .shp <- if (identical(prep$covFamily[j], "log")) "power" else "lin"
        .r <- .vaeShapeBeta(.shp, .ctr, fit$beta[k, j])
      }
      ## "." separates the pieces, matching the rest of nlmixr2's generated and
      ## conventional parameter names (eta.cl, add.sd, prop.sd) rather than "_".
      ## The categorical branch is built from the raw covariate and its level
      ## directly instead of reusing the design column name (which glues them
      ## with "_"), so the separator is consistent across both branches.
      ## The hockey arms name themselves `hockey.low` / `hockey.hi` rather than
      ## by their internal arm shape, so the pair reads as one relationship.
      .bn <- if (identical(prep$covType[j], "continuous")) {
        paste0("beta.", thName, ".", .raw, ".", .vaeShapeCoefTag(.shp))
      } else if (!is.na(prep$covLevel[j]) && nzchar(prep$covLevel[j])) {
        paste0("beta.", thName, ".", prep$covRaw[j], ".", prep$covLevel[j])
      } else {
        paste0("beta.", thName, ".", covNames[j])
      }
      ## names are built from user-facing pieces, so two distinct columns can in
      ## principle collide (a continuous WT written as "log" beside a 0/1 data
      ## column literally named WT_log); keep them distinct rather than letting
      ## one silently overwrite the other in betaVals
      .bn <- .vaeUniqueName(.bn, names(betaVals))
      icAdj[k] <- icAdj[k] + .r$interceptAdj
      .enc <- if (identical(prep$covShape[j], "cat")) {
        ## an indicator enters as written in the data (bare 0/1 column) or as an
        ## explicit level comparison -- either way it is never re-centered
        prep$covExpr[j]
      } else {
        .vaeShapeExpr(.shp, .raw, .ctr)
      }
      terms <- c(terms, paste0(.bn, " * ", .enc))
      betaVals[[.bn]] <- .r$beta
    }
    ## Inject the covariate terms FLAT (no wrapping parentheses): the mu-ref line
    ## is `p <- exp(theta + eta)`, so replacing `theta` with `theta + beta*cov`
    ## keeps the additive `exp(theta + beta*cov + eta)` form rxode2 recognizes as
    ## a mu-referenced exp() parameter.  Wrapping in parens -- `exp((theta +
    ## beta*cov) + eta)` -- hides the exp() back-transform from muRefCurEval, so
    ## the theta prints on the raw log scale instead of back-transformed.
    .termTxt <- paste(terms, collapse = " + ")
    .new <- .vaeInjectCov(.lines[[.idx]], thName, prep$etaNames[k], .termTxt)
    if (is.null(.new)) {
      ## no mu-referenced occurrence found (an unusual line shape); fall back to
      ## the FIRST textual occurrence rather than every one of them
      .new <- str2lang(sub(paste0("\\b", thName, "\\b"),
                           paste0(thName, " + ", .termTxt),
                           deparse1(.lines[[.idx]])))
    }
    ui2 <- do.call(rxode2::model, list(ui2, .new))
  }

  ## 2. set ini() estimates to the VAE solution
  .setIni <- function(u, expr) do.call(rxode2::ini, list(u, str2lang(expr)))
  for (k in seq_along(thetaNames)) {
    ## a free/fixed eta has no structural theta (thetaForEta == NA) -- its
    ## population location is already a literal in the model, so only set omega
    if (!is.na(thetaNames[k])) {
      ui2 <- .setIni(ui2, paste0(thetaNames[k], " <- ",
                                 signif(fit$zPop[k] + icAdj[k], 12)))
    }
  }
  ui2 <- .omegaWriteIni(ui2, .omegaFitMat(fit, fit$prep$etaNames))
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
  ## stale: an injected covariate-coefficient theta (beta.<par>.<cov>) is added to
  ## the iniDf as a theta but ALSO stays listed as a covariate.  The augmented
  ## covariance solve then declares that beta. both as its THETA[k] and as a
  ## phantom data covariate, and fails ("required for solving: beta....").
  ## Re-parsing the accumulated model function yields a consistent theta/covariate
  ## classification (verified: only the true data covariates remain).
  rxode2::assertRxUi(ui2$fun)
}

#' Update a PINNED VAE fit's model with the estimates.
#'
#' Unlike `.vaeUpdateModel` (which injects fresh `beta.<par>.<cov>` terms into a
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
    ## the pair pins to the column matching the form it was written in, which is
    ## how the coefficient transfers with no re-parameterization
    .j <- .vaePinColumn(prep, .inRows[.r, , drop = FALSE])
    .sel <- !is.null(fit$selected) && !is.na(.j) && isTRUE(fit$selected[.k, .j])
    .betaVal <- if (.sel) fit$beta[.k, .j] else 0
    ui2 <- .setIni(ui2, paste0(.inRows$coefName[.r], " <- ", signif(.betaVal, 12)))
  }

  ## 2. structural population thetas + omega (block-aware)
  for (k in seq_along(thetaNames)) {
    if (!is.na(thetaNames[k])) {
      ui2 <- .setIni(ui2, paste0(thetaNames[k], " <- ", signif(fit$zPop[k], 12)))
    }
  }
  ui2 <- .omegaWriteIni(ui2, .omegaFitMat(fit, prep$etaNames))

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
                   covNames = fit$covNames, zPop = fit$zPop, omega = fit$omega,
                   omegaMat = fit$omegaMat, a = fit$a,
                   covSelectMethodUsed = fit$covSelectMethodUsed,
                   seed = .control$seed)
  ## the VAE optimization walk (standard parHistData -> $parHist accessor)
  if (!is.null(fit$parHist)) .ret$parHistData <- fit$parHist
  nmObjHandleControlObject(.control, .ret) # stores $vaeControl for nmObjGetControl.vae
  .vaeControlToFoceiControl(.ret)
  ## ---- foreign-method output contract -------------------------------------
  ## `.ret` is a fresh env, so it carries none of the state the output builder
  ## needs.  Supply it the way an out-of-tree method must (the reference is
  ## babelmixr2's R/saemix.R): first the data-derived state, then the fit items.
  ## Leaving these to be derived works for a plain model but breaks on IOV -- the
  ## derivation looks a theta up per eta via muRefTable, and an IOV eta
  ## (rx.iov.<v>.<occ>) has NO muRefTable row, so the lookup is zero-length
  ## ("invalid second argument of length 0").  Every lookup below is guarded on
  ## length()==1 so occasion etas fall through harmlessly.
  if (!exists("dataSav", envir = .ret, inherits = FALSE)) {
    .foceiPreProcessData(env$data, .ret, .ui2, .ret$control$rxControl)
  }
  .idf2 <- .ui2$iniDf
  .etaU <- .idf2$name[!is.na(.idf2$neta1) & .idf2$neta1 == .idf2$neta2]
  ## 1. fullTheta -- every non-eta ini() entry, already carrying the VAE estimates
  ##    (.vaeUpdateModel wrote them into .ui2)
  if (!exists("fullTheta", envir = .ret, inherits = FALSE)) {
    .ret$fullTheta <- setNames(.idf2$est[is.na(.idf2$neta1)], .idf2$name[is.na(.idf2$neta1)])
  }
  ## 2. etaObf -- ID + one column per UI eta + OBJI, in eta order
  if (!exists("etaObf", envir = .ret, inherits = FALSE)) {
    .ids <- unique(.ret$dataSav$ID)
    .em <- .etaMat[, intersect(.etaU, colnames(.etaMat)), drop = FALSE]
    if (nrow(.em) == length(.ids)) {
      .eo <- as.data.frame(.em)
      .eo$ID <- .ids
      .eo <- .eo[, c("ID", colnames(.em)), drop = FALSE]
      .eo$OBJI <- NA_real_
      .ret$etaObf <- .eo
    }
  }
  ## 3. omega -- dimnamed by the UI eta names, from the updated iniDf: the VAE
  ##    estimates the full modeled block (diagonal + declared off-diagonals);
  ##    an occasion eta keeps whatever the model fixed it at
  if (!exists("omega", envir = .ret, inherits = FALSE)) {
    .om <- .omegaBlockFromIniDf(.idf2, .etaU)$mat
    .om[!is.finite(.om)] <- 0
    .ret$omega <- .om
  }
  ## 4/5. cov + objective are deliberately NOT set: the builder derives them from
  ##      the FOCEi inner pass at the VAE estimates, which is how the VAE reports
  ##      its objective today.  Setting them here would change every fit.
  ## 6. remaining metadata ($method/$extra/$est set above)
  if (!exists("message", envir = .ret, inherits = FALSE)) .ret$message <- ""
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
