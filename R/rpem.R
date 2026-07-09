# RPEM estimation method (design/rpem/).  This file wires the validated C++ K=1
# engine (src/rpem.cpp) to a mu-referenced rxode2 UI: it classifies the model
# parameters and drives the E-step -> conjugate M-step iteration loop.
#
# M1 scope: single endpoint, diagonal omega, one additive residual (add.sd).
# Structural parameters without a between-subject eta are held fixed for now
# (their numeric M-step update is a follow-up); mixtures/IOV/censoring are later
# milestones.

#' Classify a mu-referenced UI into the pieces the RPEM engine needs.
#'
#' @param ui rxode2 UI object.
#' @return list with the engine inputs (see body).
#' @noRd
#' Extract the TBS transform code (yj) and bounds (low/hi) from the RPEM model.
#'
#' The generated `rpemModel0` sets `rx_yj_ ~ <code>`, `rx_low_ ~ <lo>`,
#' `rx_hi_ ~ <hi>` -- the same yj0 code the model's rxTBS uses, so passing it to
#' C++ `_powerD`/`_powerDD` reproduces the model transform exactly.
#' @noRd
.rpemExtractTBS <- function(ui) {
  .m <- ui$rpemModel0
  .stmts <- as.list(.m[[2]])
  .grab <- function(nm) {
    for (.s in .stmts) {
      if (is.call(.s) && identical(.s[[1]], as.name("~")) &&
          is.name(.s[[2]]) && as.character(.s[[2]]) == nm && is.numeric(.s[[3]])) {
        return(as.numeric(.s[[3]]))
      }
    }
    NA_real_
  }
  list(yj = as.integer(.grab("rx_yj_")), low = .grab("rx_low_"), hi = .grab("rx_hi_"))
}

#' Per-observation endpoint index (0-based) in the E-step solve order.
#'
#' Mirrors SAEM's `ix_endpnt`: each observation is mapped to its endpoint via the
#' data `dvid` (preferred) or `cmt` tag, in the (id, time) order rxode2 solves.
#' @noRd
.rpemEndptIndex <- function(data, cl) {
  .nm <- names(data)
  .idCol <- .nm[tolower(.nm) == "id"][1]
  .tCol <- .nm[tolower(.nm) == "time"][1]
  .evidCol <- .nm[tolower(.nm) == "evid"][1]
  .obs <- if (!is.na(.evidCol)) data[data[[.evidCol]] == 0, , drop = FALSE] else data
  .obs <- .obs[order(.obs[[.idCol]], .obs[[.tCol]]), , drop = FALSE]
  .dvidCol <- .nm[tolower(.nm) == "dvid"][1]
  .cmtCol <- .nm[tolower(.nm) == "cmt"][1]
  if (!is.na(.dvidCol)) {
    .idx <- match(.obs[[.dvidCol]], cl$endpt$dvid)
  } else if (!is.na(.cmtCol)) {
    .idx <- match(.obs[[.cmtCol]], cl$endpt$cmt)
  } else {
    stop("RPEM multiple-endpoint data needs a 'dvid' or 'cmt' column to tag endpoints")
  }
  if (anyNA(.idx)) stop("could not map every observation to a model endpoint")
  as.integer(.idx - 1L)
}

.rpemClassify <- function(ui) {
  .ini <- ui$iniDf
  .thetas <- .ini[!is.na(.ini$ntheta), , drop = FALSE]
  .thetas <- .thetas[order(.thetas$ntheta), , drop = FALSE]
  nTheta <- nrow(.thetas)
  .etas <- .ini[!is.na(.ini$neta1) & .ini$neta1 == .ini$neta2, , drop = FALSE]
  .etas <- .etas[order(.etas$neta1), , drop = FALSE]
  nEta <- nrow(.etas)
  if (nEta == 0L) stop("RPEM requires at least one between-subject random effect")
  if (any(.ini$fix)) stop("RPEM does not yet support fixed (fix()) parameters")
  # param vector order (rpemParams): THETA[1..nTheta], ETA[1..nEta] (+ DV from data)
  base <- c(.thetas$est, rep(0, nEta))
  etaIdx <- as.integer(nTheta + seq_len(nEta) - 1L)          # 0-based
  omega0 <- diag(.etas$est, nEta)
  # mu-referenced typical-value theta for each eta (in eta order)
  .mu <- ui$muRefDataFrame
  .muName <- .mu$theta[match(.etas$name, .mu$eta)]
  if (anyNA(.muName)) stop("RPEM requires every random effect to be mu-referenced")
  muIdx <- as.integer(match(.muName, .thetas$name) - 1L)     # 0-based theta positions
  # residual: additive (add), proportional (prop), or combined (add + prop).
  .res <- .thetas[!is.na(.thetas$err), , drop = FALSE]
  .errs <- sort(.res$err)
  .pred <- ui$predDf
  .nEndpt <- nrow(.pred)
  .endpt <- NULL
  if (.nEndpt > 1L) {
    # multiple endpoints (mirror SAEM): one residual update per endpoint.  This
    # increment supports an additive or proportional residual per endpoint; the
    # E-step already computes the joint multi-endpoint likelihood.
    .iniErr <- .ini[!is.na(.ini$err), , drop = FALSE]
    .cmt <- integer(.nEndpt); .dvid <- integer(.nEndpt); .et <- integer(.nEndpt)
    .sclIdx <- integer(.nEndpt); .scl0 <- numeric(.nEndpt); .enm <- character(.nEndpt)
    for (.b in seq_len(.nEndpt)) {
      .cond <- as.character(.pred$cond[.b])
      if (!(as.character(.pred$errType[.b]) %in% c("add", "prop")))
        stop("RPEM multiple-endpoint currently supports an additive or proportional residual per endpoint")
      .er <- .iniErr[as.character(.iniErr$condition) == .cond &
                       .iniErr$err %in% c("add", "prop"), , drop = FALSE]
      if (nrow(.er) != 1L)
        stop("RPEM multiple-endpoint expects exactly one add/prop residual per endpoint")
      .cmt[.b] <- as.integer(.pred$cmt[.b])
      .dvid[.b] <- as.integer(.pred$dvid[.b])
      .et[.b] <- if (.er$err[1] == "prop") 1L else 0L
      .sclIdx[.b] <- as.integer(match(.er$name, .thetas$name) - 1L)
      .scl0[.b] <- .er$est; .enm[.b] <- .cond
    }
    errType <- 5L; errName <- "multiEndpoint"
    addSdIdx <- NA_integer_; addSd0 <- NA_real_; propSdIdx <- NA_integer_; propSd0 <- NA_real_
    .endpt <- list(nEndpt = .nEndpt, cmt = .cmt, dvid = .dvid, errType = .et,
                   sclIdx = .sclIdx, scl0 = .scl0, name = .enm)
  } else if (nrow(.res) == 1L && .res$err[1] %in% c("add", "prop")) {
    errType <- if (.res$err[1] == "prop") 1L else 0L
    # addSdIdx points at the single residual (holds add.sd or prop.sd)
    addSdIdx <- as.integer(match(.res$name, .thetas$name) - 1L)
    propSdIdx <- NA_integer_; propSd0 <- NA_real_
    errName <- .res$err[1]
    addSd0 <- .res$est
  } else if (nrow(.res) == 2L && identical(.errs, c("add", "prop"))) {
    errType <- 2L; errName <- "add+prop"
    .addRow <- .res[.res$err == "add", , drop = FALSE]
    .propRow <- .res[.res$err == "prop", , drop = FALSE]
    addSdIdx <- as.integer(match(.addRow$name, .thetas$name) - 1L)
    propSdIdx <- as.integer(match(.propRow$name, .thetas$name) - 1L)
    addSd0 <- .addRow$est; propSd0 <- .propRow$est
  } else if (any(.res$err %in% c("boxCox", "yeoJohnson"))) {
    # transform-both-sides: additive residual on the transformed scale + a
    # dynamic (estimated) Box-Cox / Yeo-Johnson lambda (errType 3).
    .lamRow <- .res[.res$err %in% c("boxCox", "yeoJohnson"), , drop = FALSE]
    .sclRow <- .res[!(.res$err %in% c("boxCox", "yeoJohnson")), , drop = FALSE]
    if (nrow(.lamRow) != 1L || nrow(.sclRow) != 1L || .sclRow$err[1] != "add")
      stop("RPEM TBS currently supports a single additive residual with a boxCox/yeoJohnson transform")
    errType <- 3L; errName <- paste0("add+", .lamRow$err[1])
    addSdIdx <- as.integer(match(.sclRow$name, .thetas$name) - 1L)
    addSd0 <- .sclRow$est
    propSdIdx <- NA_integer_; propSd0 <- NA_real_
    lambdaIdx <- as.integer(match(.lamRow$name, .thetas$name) - 1L)
    lambda0 <- .lamRow$est
    .tbs <- .rpemExtractTBS(ui)
  } else if (any(.res$err %in% c("pow", "pow2"))) {
    # power error: variance (prop.sd * cp^power)^2, both estimated (errType 4).
    .sclRow <- .res[.res$err == "pow", , drop = FALSE]
    .powRow <- .res[.res$err == "pow2", , drop = FALSE]
    if (nrow(.sclRow) != 1L || nrow(.powRow) != 1L)
      stop("RPEM power error currently supports a single pow(scale, exponent) term")
    errType <- 4L; errName <- "pow"
    addSdIdx <- as.integer(match(.sclRow$name, .thetas$name) - 1L)  # holds the scale (prop.sd)
    addSd0 <- .sclRow$est
    powIdx <- as.integer(match(.powRow$name, .thetas$name) - 1L)
    pow0 <- .powRow$est
    propSdIdx <- NA_integer_; propSd0 <- NA_real_
  } else {
    stop("RPEM currently supports additive, proportional, combined (add + prop), TBS (add + boxCox/yeoJohnson), or power residual error")
  }
  if (errType != 3L) { lambdaIdx <- NA_integer_; lambda0 <- NA_real_
    .tbs <- list(yj = NA_integer_, low = NA_real_, hi = NA_real_) }
  if (errType != 4L) { powIdx <- NA_integer_; pow0 <- NA_real_ }
  # mu2 covariates on the mu-referenced (eta) params (D22): the covariate
  # coefficients are estimated via the regression M-step, not held.
  .covDf <- ui$muRefCovariateDataFrame
  if (is.null(.covDf) || nrow(.covDf) == 0L) {
    covCoefNames <- character(0); covNames <- character(0)
  } else {
    .covDf <- .covDf[.covDf$theta %in% .muName, , drop = FALSE]
    covCoefNames <- as.character(.covDf$covariateParameter)
    covNames <- as.character(.covDf$covariate)
  }
  covCoefIdx <- as.integer(match(covCoefNames, .thetas$name) - 1L)
  list(base = base, nTheta = nTheta, nEta = nEta, etaIdx = etaIdx, omega0 = omega0,
       muIdx = muIdx, mu0 = .thetas$est[muIdx + 1L],
       addSdIdx = addSdIdx, addSd0 = addSd0, errType = errType,
       propSdIdx = propSdIdx, propSd0 = propSd0, errName = errName,
       lambdaIdx = lambdaIdx, lambda0 = lambda0,
       tbsYj = .tbs$yj, tbsLow = .tbs$low, tbsHi = .tbs$hi,
       powIdx = powIdx, pow0 = pow0, endpt = .endpt,
       covCoefNames = covCoefNames, covNames = covNames, covCoefIdx = covCoefIdx,
       covCoef0 = if (length(covCoefIdx)) .thetas$est[covCoefIdx + 1L] else numeric(0),
       thetaNames = .thetas$name, etaNames = .etas$name, muNames = .muName)
}

#' Fit a mu-referenced model with RPEM (K=1 core).
#'
#' @param ui rxode2 UI object.
#' @param data Data frame (with a DV column).
#' @param control `rpemControl()`.
#' @return list of estimates (`mu`, `omega`, `addSd`) plus per-iteration traces.
#' @noRd
.rpemFit <- function(ui, data, control = rpemControl()) {
  .cl <- .rpemClassify(ui)
  .m <- ui$rpemRxModel$predOnly
  .nm <- c(paste0("THETA[", seq_len(.cl$nTheta), "]"),
           paste0("ETA[", seq_len(.cl$nEta), "]"))
  .e <- new.env()
  .e$predOnly <- .m
  .e$rxControl <- rxode2::rxControl(atol = control$atol, rtol = control$rtol)
  .e$param <- stats::setNames(.cl$base, .nm)
  .e$data <- data

  # mu2 covariate design (D22): for a single random effect, estimate the typical
  # value + covariate coefficients via the regression M-step.  design row i =
  # [1, cov1_i, ...] with per-subject covariate values (solve order = sorted id).
  .useReg <- (.cl$nEta == 1L)
  if (!.useReg && length(.cl$covCoefNames) > 0L)
    stop("RPEM does not yet support covariates with more than one random effect")
  if (.useReg) {
    .idCol <- names(data)[tolower(names(data)) == "id"][1]
    if (is.na(.idCol)) stop("data must have an ID column")
    .ids <- sort(unique(data[[.idCol]]))
    .design <- matrix(1, length(.ids), 1L + length(.cl$covNames))
    for (.j in seq_along(.cl$covNames)) {
      .cv <- .cl$covNames[.j]
      .design[, .j + 1L] <- vapply(.ids, function(i) data[[.cv]][data[[.idCol]] == i][1], numeric(1))
    }
    coefs <- c(.cl$mu0, .cl$covCoef0)
  }

  .comb <- (.cl$errType == 2L)
  .tbs <- (.cl$errType == 3L)
  .pow <- (.cl$errType == 4L)
  .multi <- (.cl$errType == 5L)
  base <- .cl$base; mu <- .cl$mu0; omega <- .cl$omega0
  addSd <- .cl$addSd0; propSd <- if (.comb) .cl$propSd0 else NA_real_
  lambda <- if (.tbs) .cl$lambda0 else NA_real_
  power <- if (.pow) .cl$pow0 else NA_real_
  if (.multi) { .endptIdx <- .rpemEndptIndex(data, .cl); sdVec <- .cl$endpt$scl0 }
  niter <- control$niter
  coefTr <- if (.useReg) matrix(0, niter, ncol(.design)) else NULL
  muTr <- matrix(0, niter, .cl$nEta); omTr <- matrix(0, niter, .cl$nEta)
  sdTr <- numeric(niter); propTr <- numeric(niter); lamTr <- numeric(niter)
  powTr <- numeric(niter); llTr <- numeric(niter)
  sdMat <- if (.multi) matrix(0, niter, .cl$endpt$nEndpt) else NULL
  for (.it in seq_len(niter)) {
    if (.useReg) {
      base[.cl$muIdx + 1L] <- coefs[1]
      if (length(.cl$covCoefIdx)) base[.cl$covCoefIdx + 1L] <- coefs[-1]
    } else {
      base[.cl$muIdx + 1L] <- mu
    }
    if (!.multi) base[.cl$addSdIdx + 1L] <- addSd
    if (.comb) base[.cl$propSdIdx + 1L] <- propSd
    if (.tbs) base[.cl$lambdaIdx + 1L] <- lambda
    if (.pow) base[.cl$powIdx + 1L] <- power
    if (.multi) base[.cl$endpt$sclIdx + 1L] <- sdVec
    rxode2::rxSetSeed(control$seed + .it)
    .est <- rpemEstepK1Draw(.e, base, .cl$etaIdx, omega, control$nGauss, control$cores)
    if (.comb) {
      .ms <- rpemMstepK1Comb(.design, coefs, addSd, propSd, control$nMH, control$mhBurn)
      coefs <- .ms$coefs; omega <- matrix(.ms$omega, 1, 1)
      addSd <- .ms$addSd; propSd <- .ms$propSd
      coefTr[.it, ] <- coefs; muTr[.it, ] <- coefs[1]; omTr[.it, ] <- .ms$omega
    } else if (.tbs) {
      .ms <- rpemMstepK1TBS(.design, coefs, addSd, lambda, .cl$tbsYj,
                            .cl$tbsLow, .cl$tbsHi, control$nMH, control$mhBurn)
      coefs <- .ms$coefs; omega <- matrix(.ms$omega, 1, 1)
      addSd <- .ms$addSd; lambda <- .ms$lambda
      coefTr[.it, ] <- coefs; muTr[.it, ] <- coefs[1]; omTr[.it, ] <- .ms$omega
    } else if (.pow) {
      .ms <- rpemMstepK1Pow(.design, coefs, addSd, power, control$nMH, control$mhBurn)
      coefs <- .ms$coefs; omega <- matrix(.ms$omega, 1, 1)
      addSd <- .ms$propSd; power <- .ms$power
      coefTr[.it, ] <- coefs; muTr[.it, ] <- coefs[1]; omTr[.it, ] <- .ms$omega
    } else if (.multi) {
      .ms <- rpemMstepK1Multi(.design, coefs, .endptIdx, .cl$endpt$errType,
                              control$nMH, control$mhBurn)
      coefs <- .ms$coefs; omega <- matrix(.ms$omega, 1, 1); sdVec <- .ms$sd
      coefTr[.it, ] <- coefs; muTr[.it, ] <- coefs[1]; omTr[.it, ] <- .ms$omega
      sdMat[.it, ] <- sdVec
    } else if (.useReg) {
      .ms <- rpemMstepK1Reg(.design, coefs, .cl$errType, control$nMH, control$mhBurn)
      coefs <- .ms$coefs; omega <- matrix(.ms$omega, 1, 1); addSd <- .ms$addSd
      coefTr[.it, ] <- coefs; muTr[.it, ] <- coefs[1]; omTr[.it, ] <- .ms$omega
    } else {
      .ms <- rpemMstepK1(mu, addSd, control$nMH, control$mhBurn)
      mu <- .ms$mu; omega <- .ms$omega; addSd <- .ms$addSd
      muTr[.it, ] <- mu; omTr[.it, ] <- diag(omega)
    }
    sdTr[.it] <- addSd; propTr[.it] <- if (.comb) propSd else NA_real_
    lamTr[.it] <- if (.tbs) lambda else NA_real_
    powTr[.it] <- if (.pow) power else NA_real_
    llTr[.it] <- .est$lnL
  }

  # Final estimate = mean over the converged iterations.
  .k <- min(control$collect, niter)
  .w <- (niter - .k + 1L):niter
  muHat <- colMeans(muTr[.w, , drop = FALSE])
  omHat <- colMeans(omTr[.w, , drop = FALSE])
  sdHat <- mean(sdTr[.w])
  propHat <- if (.comb) mean(propTr[.w]) else NA_real_
  lambdaHat <- if (.tbs) mean(lamTr[.w]) else NA_real_
  powerHat <- if (.pow) mean(powTr[.w]) else NA_real_
  endptSdHat <- if (.multi) colMeans(sdMat[.w, , drop = FALSE]) else NULL
  covCoefHat <- if (.useReg && length(.cl$covCoefIdx))
                  colMeans(coefTr[.w, -1, drop = FALSE]) else numeric(0)

  # One final E-step at the converged estimates to compute per-subject EBEs
  # (posterior-mean etas, Eq 53): EBE_i = sum_j eta_ij * w_ij with the
  # self-normalized importance weights w_ij = softmax_j(log p_ij).
  base[.cl$muIdx + 1L] <- muHat
  if (length(covCoefHat)) base[.cl$covCoefIdx + 1L] <- covCoefHat
  if (!.multi) base[.cl$addSdIdx + 1L] <- sdHat
  if (.comb) base[.cl$propSdIdx + 1L] <- propHat
  if (.tbs) base[.cl$lambdaIdx + 1L] <- lambdaHat
  if (.pow) base[.cl$powIdx + 1L] <- powerHat
  if (.multi) base[.cl$endpt$sclIdx + 1L] <- endptSdHat
  omegaHat <- diag(omHat, .cl$nEta)
  rxode2::rxSetSeed(control$seed)
  .fe <- rpemEstepK1Draw(.e, base, .cl$etaIdx, omegaHat, control$nGauss, control$cores)
  rpemFree()
  .nG <- control$nGauss
  .nsub <- length(.fe$logp) / .nG
  .etaM <- matrix(.fe$eta, ncol = .cl$nEta)
  ebe <- matrix(0, .nsub, .cl$nEta, dimnames = list(NULL, .cl$etaNames))
  for (.i in seq_len(.nsub)) {
    .idx <- ((.i - 1L) * .nG + 1L):(.i * .nG)
    .lw <- .fe$logp[.idx]; .wt <- exp(.lw - max(.lw)); .wt <- .wt / sum(.wt)
    ebe[.i, ] <- colSums(.etaM[.idx, , drop = FALSE] * .wt)
  }

  list(mu = stats::setNames(muHat, .cl$muNames),
       omega = stats::setNames(omHat, .cl$etaNames),
       addSd = sdHat, propSd = propHat, lambda = lambdaHat, power = powerHat,
       endptSd = endptSdHat,
       covCoef = stats::setNames(covCoefHat, .cl$covCoefNames),
       ebe = ebe,
       lnL = llTr, muTrace = muTr, omegaTrace = omTr, sdTrace = sdTr,
       classify = .cl)
}

#' Validate the RPEM control (est="rpem")
#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.rpem <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- rpemControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("rpemControl", .ctl)
  if (!inherits(.ctl, "rpemControl")) {
    .minfo("invalid control for `est=\"rpem\"`, using default")
    .ctl <- rpemControl()
  } else {
    .ctl <- do.call(rpemControl, .ctl)
  }
  .ctl
}

#' Assemble a full nlmixr2FitData from RPEM estimates via eval-only FOCEI.
#'
#' Mirrors SAEM's finalize: set the estimated theta/omega and the RPEM EBEs on
#' the UI, then let FOCEI evaluate (0 outer/inner iterations) at those fixed
#' values to compute EBEs/residuals/tables (design/rpem/09).
#' @noRd
.rpemBuildFit <- function(env, ui, control, rfit) {
  .ret <- new.env(parent = emptyenv())
  .ret$table <- env$table
  .rxControl <- rxode2::rxControl(atol = control$atol, rtol = control$rtol,
                                  method = "liblsoda")
  .foceiPreProcessData(env$data, .ret, ui, .rxControl)
  .cl <- rfit$classify
  # M1 holds non-mu-ref structural params (no eta, not the residual) at their ini
  # values; mark them fixed so the fit reports them as held rather than assigning
  # FOCEI-covariance SEs.  They get a proper numeric M-step update later (D20).
  .resNames <- if (.cl$errType == 5L) .cl$thetaNames[.cl$endpt$sclIdx + 1L]
               else .cl$thetaNames[.cl$addSdIdx + 1L]
  if (.cl$errType == 2L) .resNames <- c(.resNames, .cl$thetaNames[.cl$propSdIdx + 1L])
  if (.cl$errType == 3L) .resNames <- c(.resNames, .cl$thetaNames[.cl$lambdaIdx + 1L])
  if (.cl$errType == 4L) .resNames <- c(.resNames, .cl$thetaNames[.cl$powIdx + 1L])
  .heldNames <- setdiff(.cl$thetaNames,
                        c(.cl$muNames, .resNames, .cl$covCoefNames))
  if (length(.heldNames) > 0L) {
    .uiD <- rxode2::rxUiDecompress(ui)
    .idf <- .uiD$iniDf
    .idf$fix[.idf$name %in% .heldNames] <- TRUE
    assign("iniDf", .idf, envir = .uiD)
    ui <- rxode2::rxUiCompress(.uiD)
  }
  .ret$ui <- ui
  # full theta (named over all theta names): mu-referenced -> RPEM mu, additive
  # residual -> RPEM add.sd, held structural -> ini values.
  .tn <- .cl$thetaNames
  .ft <- stats::setNames(.cl$base[seq_along(.tn)], .tn)
  .ft[.cl$muNames] <- rfit$mu
  if (.cl$errType == 5L) .ft[.tn[.cl$endpt$sclIdx + 1L]] <- rfit$endptSd
  else .ft[.tn[.cl$addSdIdx + 1L]] <- rfit$addSd
  if (.cl$errType == 2L) .ft[.tn[.cl$propSdIdx + 1L]] <- rfit$propSd
  if (.cl$errType == 3L) .ft[.tn[.cl$lambdaIdx + 1L]] <- rfit$lambda
  # power: exponent set here; the scale (prop.sd) is set via the addSd slot above.
  if (.cl$errType == 4L) .ft[.tn[.cl$powIdx + 1L]] <- rfit$power
  if (length(rfit$covCoef)) .ft[.cl$covCoefNames] <- rfit$covCoef
  .ret$fullTheta <- .ft
  # omega (diagonal) named over etas
  .om <- diag(rfit$omega, .cl$nEta)
  dimnames(.om) <- list(.cl$etaNames, .cl$etaNames)
  .ret$omega <- .om
  # per-subject EBEs -> etaMat for the FOCEI eval
  .eb <- rfit$ebe
  colnames(.eb) <- .cl$etaNames
  .ret$.etaMat <- .eb
  .ret$.etaMatBase <- .eb
  .ret$etaObf <- data.frame(ID = seq_len(nrow(.eb)),
                            stats::setNames(as.data.frame(.eb), .cl$etaNames),
                            OBJI = NA)
  .nlmixr2FitUpdateParams(.ret)
  # covMethod "r,s" makes FOCEI compute the covariance (hence SEs / %RSE / CIs)
  # at the fixed RPEM estimates.  These are FOCEI-covariance-based SEs; the
  # paper's Fisher-score SEs (design/rpem/08) are a later refinement.
  .foceiControl <- foceiControl(maxOuterIterations = 0L, maxInnerIterations = 0L,
                                covMethod = "r,s", etaMat = .eb, scaleTo = 0,
                                calcTables = TRUE, interaction = 1L,
                                rxControl = .rxControl, est = "rpem")
  .ret$control <- .foceiControl
  .ret$est <- "rpem"
  .ret$ofvType <- "rpem"
  .ret$adjObf <- TRUE
  # Store the control on the fit env so nmObjGetControl (hence methodOde, etc.)
  # resolve later.  nmObjHandleControlObject may null env$control, so pass the
  # saved .foceiControl variable (not .ret$control) to nlmixr2CreateOutputFromUi.
  nmObjHandleControlObject(.foceiControl, .ret)
  # Provide the compiled FOCEI model so the residual-table getters resolve:
  # nmObjGet.innerModel -> foceiModel$inner, nmObjGetIpredModel.default ->
  # foceiModel$predOnly (est="rpem" has no dedicated ipred-model getter). dataSav
  # is recomputed from origData on demand.
  .ret$foceiModel <- ui$focei
  .out <- nlmixr2CreateOutputFromUi(.ret$ui, data = .ret$origData,
                                    control = .foceiControl, table = .ret$table,
                                    env = .ret, est = "rpem")
  # The eval-only FOCEI pass produces the full nlmixr2FitData (parFixedDf, objDf,
  # omega, EBEs, shrinkage, and the per-observation IPRED/PRED/CWRES table).  If
  # anything goes wrong, require the full data.frame and let nlmixr2Est.rpem fall
  # back to the lightweight estimates object.
  if (!inherits(.out, "nlmixr2FitData")) {
    stop("rpem residual-table step incomplete")
  }
  .out
}

#' Retrieve the (focei eval) control stored on an RPEM fit.
#' @rdname nmObjGetControl
#' @export
nmObjGetControl.rpem <- function(x, ...) {
  .env <- x[[1]]
  for (.name in c("foceiControl0", "control")) {
    if (exists(.name, .env)) {
      .control <- get(.name, .env)
      if (inherits(.control, "foceiControl")) return(.control)
    }
  }
  stop("cannot find rpem control object", call.=FALSE)
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.rpem <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (!inherits(.control, "rpemControl")) .control <- rpemControl()
  # M1: single-endpoint, mu-referenced, diagonal omega, additive residual.
  # .rpemClassify raises a clear error if the model is outside that scope.
  .fit <- .rpemFit(.ui, env$data, .control)
  # Try the full nlmixr2FitData via eval-only FOCEI; fall back to the lightweight
  # estimates object if that path errors, so est="rpem" always returns something.
  .full <- tryCatch(.rpemBuildFit(env, .ui, .control, .fit),
                    error = function(e) {
                      .minfo(paste0("rpem: full fit object unavailable (",
                                    conditionMessage(e), "); returning estimates"))
                      NULL
                    })
  if (!is.null(.full)) return(.full)
  .fit$ui <- .ui
  .fit$control <- .control
  class(.fit) <- "nlmixr2rpem"
  .fit
}

#' @export
print.nlmixr2rpem <- function(x, ...) {
  cat("RPEM fit (K=1 minimal core -- not yet a full nlmixr2 fit object)\n")
  cat("-- typical values (mu):\n"); print(round(x$mu, 4))
  cat("-- omega (between-subject variance):\n"); print(round(x$omega, 4))
  cat(sprintf("-- add.sd: %.4f\n", x$addSd))
  cat(sprintf("-- final lnL: %.3f over %d iterations\n",
              x$lnL[length(x$lnL)], length(x$lnL)))
  invisible(x)
}
