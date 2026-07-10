# vaeOutput.R -- build the fitted nlmixr2 output from a trained VAE. The model is
# updated with the selected covariate effects (exact centered expressions, e.g.
# ka <- exp(lka + beta_lka_WT*log(WT/79.6) + eta.ka)) and the ini() estimates are
# set to the VAE solution. The objective (FOCE-linearized marginal -2LL with
# M2/M3/M4 censoring and MAP-EBE re-optimization) and the linearization-Hessian
# covariance are computed by the C++/OpenMP kernel vaeFoceLik; then
# nlmixr2CreateOutputFromUi assembles the standard nlmixr2FitData (parFixed/SEs,
# EBEs, residuals, tables) WITHOUT running focei estimation (a 0-iteration focei
# control drives only the residual/table engine at the supplied etas).

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
  if (nrow(.errRow) > 0L)
    ui2 <- .setIni(ui2, paste0(.errRow$name[1], " <- ", signif(fit$a, 12)))
  ui2
}

#' Flatten the decoder linearization for the C++ likelihood: one solve per
#' subject at the encoder z gives f0, J0 = d f/d eta and the residual variance R
#' for every observation, stacked over subjects (subjOffset marks the blocks).
#' zLin = encoder posterior means (the linearization point / individual params).
#' @noRd
.vaeLinPrecomp <- function(fit) {
  prep <- fit$prep; N <- prep$N; zDim <- prep$zDim
  baseline <- fit$zPop
  th <- .vaeBuildTh(prep, baseline, fit$a)
  f0 <- numeric(0); Rv <- numeric(0); y <- numeric(0)
  cens <- integer(0); limit <- numeric(0)
  Jl <- vector("list", N); off <- integer(N + 1L)
  for (i in seq_len(N)) {
    s <- prep$subj[[i]]
    E <- .vaeDecoderSolveSubject(fit$am, th, fit$mu[i, ] - baseline, s$ev, s$times)
    if (is.null(E)) return(NULL)
    off[i + 1L] <- off[i] + s$n
    f0 <- c(f0, E$f); Rv <- c(Rv, E$R); y <- c(y, s$y)
    cens <- c(cens, s$cens); limit <- c(limit, s$limit)
    Jl[[i]] <- E$a
  }
  list(f0 = f0, J0 = do.call(rbind, Jl), R = Rv, y = y, cens = cens, limit = limit,
       subjOffset = off, zLin = fit$mu, N = N, zDim = zDim)
}

#' FOCE-linearized marginal -2LL (with M2/M3/M4 censoring, MAP-EBE, OpenMP) via
#' the C++ kernel. `zPopMat` is the per-subject population center; `Rscale`
#' multiplies the base variance (a/aFit)^2 for additive.
#' @noRd
.vaeObjective <- function(pc, zPopMat, omega, Rscale = 1, cores = 1L,
                          nMix = 1L, mixProb = 1) {
  vaeFoceLik(pc$f0, pc$J0, pc$R, pc$y, pc$cens, pc$limit, pc$zLin, zPopMat, omega,
             as.integer(pc$subjOffset), Rscale, 0L, 30L, 1e-8, as.integer(cores),
             as.integer(nMix), as.numeric(mixProb))
}

#' Thread count for the VAE likelihood: from the rxode2 control object's `cores`,
#' falling back to the rxode2 global thread count.
#' @noRd
.vaeCores <- function(control) {
  .c <- tryCatch(control$rxControl$cores, error = function(e) NULL)
  if (is.null(.c) || is.na(.c) || .c < 1L)
    .c <- tryCatch(rxode2::getRxThreads(), error = function(e) 1L)
  as.integer(.c)
}

#' Linearization-Hessian covariance of the fixed effects (intercepts + covariate
#' betas + residual a). The -2LL is an explicit function of these given the fixed
#' decoder f/J (no ODE re-solves), so the observed information is a cheap
#' numerical Hessian: cov = 2 * solve(Hessian of -2LL), holding omega.
#' Returns list(cov, names) with names matching the free-theta order, or NULL.
#' @noRd
.vaeCov <- function(fit, ui2, precomp, cores = 1L) {
  prep <- fit$prep
  .idf <- ui2$iniDf
  .thR <- .idf[!is.na(.idf$ntheta) & !(!is.na(.idf$fix) & .idf$fix), , drop = FALSE]
  .thR <- .thR[order(.thR$ntheta), , drop = FALSE]
  nm <- .thR$name; np <- length(nm)
  .map <- .foceiEtaThetaMap(ui2)
  thetaForEta <- .map$thetaForEta
  ## deterministic beta-name -> (k, j) lookup from the selected mask
  betaLU <- list()
  if (!is.null(fit$selected)) for (k in seq_len(prep$zDim)) for (j in which(fit$selected[k, ]))
    betaLU[[paste0("beta_", thetaForEta[k], "_", prep$covNames[j])]] <- c(k = k, j = j)
  base <- .thR$est
  errName <- {r <- .idf[!is.na(.idf$err) & !is.na(.idf$ntheta), ]; if (nrow(r)) r$name[1] else NA}
  obj <- function(v) {
    intercept <- setNames(fit$zPop, thetaForEta)
    aVal <- fit$a
    zPopMat <- matrix(0, prep$N, prep$zDim)
    for (p in seq_len(np)) {
      nmp <- nm[p]
      if (nmp %in% thetaForEta) intercept[nmp] <- v[p]
      else if (!is.na(errName) && nmp == errName) aVal <- v[p]
    }
    for (k in seq_len(prep$zDim)) zPopMat[, k] <- intercept[thetaForEta[k]]
    for (p in seq_len(np)) {
      bc <- betaLU[[nm[p]]]
      if (!is.null(bc)) zPopMat[, bc["k"]] <- zPopMat[, bc["k"]] + v[p] * prep$covMat[, bc["j"]]
    }
    .vaeObjective(precomp, zPopMat, fit$omega, (aVal / fit$a)^2, cores)$objective
  }
  H <- matrix(0, np, np); h <- pmax(1e-4, abs(base) * 1e-4)
  f0 <- obj(base)
  for (i in seq_len(np)) for (j in i:np) {
    if (i == j) {
      vp <- base; vp[i] <- vp[i] + h[i]; vm <- base; vm[i] <- vm[i] - h[i]
      H[i, i] <- (obj(vp) - 2 * f0 + obj(vm)) / h[i]^2
    } else {
      vpp <- base; vpp[i] <- vpp[i] + h[i]; vpp[j] <- vpp[j] + h[j]
      vpm <- base; vpm[i] <- vpm[i] + h[i]; vpm[j] <- vpm[j] - h[j]
      vmp <- base; vmp[i] <- vmp[i] - h[i]; vmp[j] <- vmp[j] + h[j]
      vmm <- base; vmm[i] <- vmm[i] - h[i]; vmm[j] <- vmm[j] - h[j]
      H[i, j] <- H[j, i] <- (obj(vpp) - obj(vpm) - obj(vmp) + obj(vmm)) / (4 * h[i] * h[j])
    }
  }
  if (!all(is.finite(H))) return(NULL)
  cov <- tryCatch(2 * solve(H), error = function(e) NULL)
  ## a non-PD linearization Hessian (e.g. a not-fully-converged fit) -> no SEs
  if (is.null(cov) || !all(is.finite(diag(cov))) || any(diag(cov) < 0)) return(NULL)
  dimnames(cov) <- list(nm, nm)
  list(cov = cov, names = nm, objective = f0)
}

#' A zero-iteration focei control for the VAE output engine. nlmixr2's residual /
#' table generation (CWRES etc.) runs through the focei engine at the SUPPLIED
#' etas (maxInner=maxOuter=0, covMethod=0) -- it assembles tables, it does NOT
#' re-estimate. Mirrors .nlmerControlToFoceiControl.
#' @noRd
.vaeControlToFoceiControl <- function(env, assign = TRUE) {
  .c <- env$vaeControl
  .fc <- foceiControl(rxControl = .c$rxControl, maxOuterIterations = 0L,
                      maxInnerIterations = 0L, covMethod = 0L, etaMat = env$etaMat,
                      sumProd = .c$sumProd, optExpression = .c$optExpression,
                      literalFix = .c$literalFix, literalFixRes = FALSE, scaleTo = 0,
                      calcTables = .c$calcTables, addProp = .c$addProp,
                      interaction = 1L, compress = .c$compress, ci = .c$ci,
                      sigdigTable = .c$sigdigTable, indTolRelax = TRUE)
  if (assign) env$control <- .fc
  .fc
}

#' Assemble the standard nlmixr2FitData from a trained VAE WITHOUT running focei
#' estimation. The model is updated with the selected covariate effects and the
#' VAE estimates; the linearization -2LL and linearization-Hessian covariance are
#' computed from the VAE; then nlmixr2CreateOutputFromUi assembles the standard
#' fit (parFixed/SEs, residuals, tables) at the VAE solution, with the encoder
#' etas supplied as the EBEs. Follows babelmixr2 nlmer.R. The ORIGINAL (pre-
#' covariate) ui is stashed in $iniDf0 for the iniUi/iniDf0 accessors.
#' @noRd
.vaeToFit <- function(env, fit) {
  .ui <- env$ui
  .control <- env$vaeControl
  .ui2 <- .vaeUpdateModel(.ui, fit)
  fit$am <- if (is.null(fit$am)) .vaeDecoderModel(.ui) else fit$am
  prep <- fit$prep

  ## VAE FOCE-linearized -2LL (C++/OpenMP, M2/M3/M4 censoring, MAP-EBE) +
  ## linearization-Hessian covariance (no ODE re-solves in the Hessian)
  precomp <- .vaeLinPrecomp(fit)
  .cores <- .vaeCores(.control)
  .obj <- .vaeObjective(precomp, fit$zPopMat, fit$omega, 1, .cores)
  .cov <- .vaeCov(fit, .ui2, precomp, .cores)

  ## assemble the output env (mirrors babelmixr2 nlmer .nlmerFamilyFit)
  .ret <- new.env(parent = emptyenv())
  .ret$table <- env$table
  .foceiPreProcessData(env$data, .ret, .ui2, .control$rxControl)
  .ret$ui <- .ui2
  .ret$adjObf <- .control$adjObf
  .idf <- .ui2$iniDf
  .thR <- .idf[!is.na(.idf$ntheta), , drop = FALSE]; .thR <- .thR[order(.thR$ntheta), ]
  .ret$fullTheta <- setNames(.thR$est, .thR$name)
  if (!is.null(.cov)) { .ret$cov <- .cov$cov; .ret$covMethod <- "linear" }
  ## MAP-refined EBEs from the Laplace re-optimization (zStar), else encoder means
  .etaMat <- .obj$zStar - fit$zPopMat
  colnames(.etaMat) <- prep$etaNames
  .ret$etaMat <- .etaMat
  .ret$etaObf <- data.frame(ID = seq_len(nrow(.etaMat)), as.data.frame(.etaMat),
                            OBJI = .obj$obji)
  .ret$omega <- diag(fit$omega, prep$zDim); dimnames(.ret$omega) <- list(prep$etaNames, prep$etaNames)
  .ret$control <- .control
  .ret$extra <- " by variational autoencoder (VAE)"
  .nlmixr2FitUpdateParams(.ret)
  nmObjHandleControlObject(.ret$control, .ret)
  if (exists("control", .ui2)) rm(list = "control", envir = .ui2)
  .ret$est <- "vae"
  .ret$objective <- .obj$objective
  .ret$model <- .ui2$ebe
  .ret$ofvType <- "vae"
  .vaeControlToFoceiControl(.ret)
  .ret$theta <- .ret$ui$saemThetaDataFrame
  .fit <- nlmixr2CreateOutputFromUi(.ret$ui, data = .ret$origData, control = .ret$control,
                                    table = .ret$table, env = .ret, est = "vae")
  .env <- .fit$env
  .env$method <- "vae"
  .env$vae <- list(elboTrace = fit$elboTrace, beta = fit$beta, selected = fit$selected,
                   covNames = fit$covNames, zPop = fit$zPop, omega = fit$omega, a = fit$a)
  ## the model STRUCTURE changed (covariate selection); stash the ORIGINAL ui so
  ## $uiIni / $iniDf0 report the original model, $finalUi the covariate model
  .env$iniDf0 <- rxode2::rxUiCompress(rxode2::rxUiDecompress(.ui))
  .fit
}
