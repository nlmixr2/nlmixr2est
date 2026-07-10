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

#' Extract the mix() mixture structure (split-ETA) from the UI.
#'
#' A model like `ka <- mix(exp(tka1 + eta.ka), p1, exp(tka2 + eta.ka))` is a
#' K-component mixture whose components differ only in the typical value
#' (tka1/tka2) and share the random effect (eta.ka).  Returns the number of
#' components, the per-component typical-value theta names (the paper's mu_k), the
#' mixed eta, and the K-1 mixture-probability theta names (P of components 1..K-1).
#' The rxode2 model selects the component via setIndMixest, so RPEM forces each
#' component in the E-step; the model uses that component's typical value.
#' @noRd
.rpemMixInfo <- function(ui, thetaNames, etaNames) {
  if (length(ui$mixProbs) == 0L) return(NULL)
  .mixCalls <- do.call(c, lapply(ui$lstExpr, .findMixCalls))
  if (length(.mixCalls) != 1L)
    stop("RPEM currently supports a single mix() call (one mixed parameter)")
  .args <- as.list(.mixCalls[[1]])[-1]
  .K <- (length(.args) + 1L) / 2L
  .compExprs <- .args[seq(1, length(.args), by = 2)]
  .probExprs <- if (.K > 1L) .args[seq(2, length(.args), by = 2)] else list()
  .muNames <- character(.K); .etaName <- NULL
  for (.k in seq_len(.K)) {
    .th <- intersect(all.vars(.compExprs[[.k]]), thetaNames)
    .et <- intersect(all.vars(.compExprs[[.k]]), etaNames)
    if (length(.th) != 1L)
      stop("RPEM mix() component must reference exactly one typical-value parameter")
    .muNames[.k] <- .th
    if (length(.et) == 1L) .etaName <- .et
  }
  if (is.null(.etaName)) stop("RPEM mix() components must share one random effect")
  .probNames <- vapply(.probExprs, function(e) as.character(e), character(1))
  list(K = as.integer(.K), muNames = .muNames, etaName = .etaName, probNames = .probNames)
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
  # param vector order (rpemParams): THETA[1..nTheta], ETA[1..nEta] (+ DV from data)
  base <- c(.thetas$est, rep(0, nEta))
  etaIdx <- as.integer(nTheta + seq_len(nEta) - 1L)          # 0-based
  omega0 <- diag(.etas$est, nEta)
  # mixture (mix(), split-ETA): the mixed eta's typical value is per-component and
  # selected by the model, so it is exempt from the plain mu-reference requirement.
  .mixI <- .rpemMixInfo(ui, .thetas$name, .etas$name)
  # mu-referenced typical-value theta for each eta (in eta order)
  .mu <- ui$muRefDataFrame
  .muName <- .mu$theta[match(.etas$name, .mu$eta)]
  .mix <- NULL
  if (!is.null(.mixI)) {
    # the mixed eta uses component 1's typical value as its nominal mu (the others
    # come from setIndMixest); non-mixed etas must still be mu-referenced.
    .mixEtaRow <- which(.etas$name == .mixI$etaName)
    .muName[.mixEtaRow] <- .mixI$muNames[1]
    if (nEta != 1L)
      stop("RPEM mixtures currently support a single mixed random effect")
    .muCompIdx <- as.integer(match(.mixI$muNames, .thetas$name) - 1L)   # 0-based tka_k
    .probIdx <- as.integer(match(.mixI$probNames, .thetas$name) - 1L)   # 0-based p_k
    .p <- .thetas$est[.probIdx + 1L]
    .mix <- list(K = .mixI$K, etaRow = as.integer(.mixEtaRow), muCompIdx = .muCompIdx,
                 muComp0 = .thetas$est[.muCompIdx + 1L], probIdx = .probIdx,
                 probNames = .mixI$probNames, muNames = .mixI$muNames,
                 w0 = c(.p, 1 - sum(.p)))
  }
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
    .sclIdx <- integer(.nEndpt); .scl0 <- numeric(.nEndpt)
    .propIdx <- integer(.nEndpt); .prop0 <- numeric(.nEndpt); .enm <- character(.nEndpt)
    for (.b in seq_len(.nEndpt)) {
      .cond <- as.character(.pred$cond[.b])
      .er <- .iniErr[as.character(.iniErr$condition) == .cond, , drop = FALSE]
      .es <- sort(.er$err)
      .cmt[.b] <- as.integer(.pred$cmt[.b]); .dvid[.b] <- as.integer(.pred$dvid[.b])
      .enm[.b] <- .cond; .propIdx[.b] <- NA_integer_; .prop0[.b] <- NA_real_
      if (nrow(.er) == 1L && .er$err[1] %in% c("add", "prop")) {
        .et[.b] <- if (.er$err[1] == "prop") 1L else 0L
        .sclIdx[.b] <- as.integer(match(.er$name, .thetas$name) - 1L); .scl0[.b] <- .er$est
      } else if (nrow(.er) == 1L && .er$err[1] == "lnorm") {
        # lognormal: additive residual on the log scale (fixed transform, no lambda)
        .et[.b] <- 6L
        .sclIdx[.b] <- as.integer(match(.er$name, .thetas$name) - 1L); .scl0[.b] <- .er$est
      } else if (nrow(.er) == 2L && identical(.es, c("add", "prop"))) {
        .et[.b] <- 2L
        .aRow <- .er[.er$err == "add", , drop = FALSE]; .pRow <- .er[.er$err == "prop", , drop = FALSE]
        .sclIdx[.b] <- as.integer(match(.aRow$name, .thetas$name) - 1L); .scl0[.b] <- .aRow$est
        .propIdx[.b] <- as.integer(match(.pRow$name, .thetas$name) - 1L); .prop0[.b] <- .pRow$est
      } else if (nrow(.er) == 2L && identical(.es, c("pow", "pow2"))) {
        # power error on this endpoint: variance (scale * cp^exponent)^2; the
        # second slot (propIdx/prop0) carries the estimated exponent.
        .et[.b] <- 4L
        .sRow <- .er[.er$err == "pow", , drop = FALSE]; .pwRow <- .er[.er$err == "pow2", , drop = FALSE]
        .sclIdx[.b] <- as.integer(match(.sRow$name, .thetas$name) - 1L); .scl0[.b] <- .sRow$est
        .propIdx[.b] <- as.integer(match(.pwRow$name, .thetas$name) - 1L); .prop0[.b] <- .pwRow$est
      } else if (nrow(.er) == 2L && "add" %in% .er$err &&
                 any(.er$err %in% c("boxCox", "yeoJohnson"))) {
        # TBS on this endpoint: additive on the transformed scale + dynamic lambda
        # (errType 3); the second slot (propIdx/prop0) carries the estimated lambda.
        .et[.b] <- 3L
        .aRow <- .er[.er$err == "add", , drop = FALSE]
        .lRow <- .er[.er$err %in% c("boxCox", "yeoJohnson"), , drop = FALSE]
        .sclIdx[.b] <- as.integer(match(.aRow$name, .thetas$name) - 1L); .scl0[.b] <- .aRow$est
        .propIdx[.b] <- as.integer(match(.lRow$name, .thetas$name) - 1L); .prop0[.b] <- .lRow$est
      } else {
        stop("RPEM multiple-endpoint supports additive, proportional, lognormal, combined (add + prop), power (pow), or TBS (add + boxCox/yeoJohnson) error per endpoint")
      }
    }
    errType <- 5L; errName <- "multiEndpoint"
    addSdIdx <- NA_integer_; addSd0 <- NA_real_; propSdIdx <- NA_integer_; propSd0 <- NA_real_
    .endpt <- list(nEndpt = .nEndpt, cmt = .cmt, dvid = .dvid, errType = .et,
                   sclIdx = .sclIdx, scl0 = .scl0,
                   propIdx = .propIdx, prop0 = .prop0, name = .enm)
  } else if (nrow(.res) == 1L && .res$err[1] %in% c("add", "prop", "lnorm")) {
    # additive, proportional, or lognormal (additive on log scale, errType 6)
    errType <- switch(.res$err[1], prop = 1L, lnorm = 6L, 0L)
    # addSdIdx points at the single residual (holds add.sd, prop.sd, or lnorm sd)
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
  # non-mu-referenced structural fixed effects (the paper's beta): thetas that are
  # not mu-referenced, not a residual/transform param, not a covariate coefficient,
  # and not user-fixed.  These are estimated by the numeric re-solve M-step; fixed
  # (fix()) thetas are held.
  .resIdx <- stats::na.omit(c(addSdIdx, propSdIdx, lambdaIdx, powIdx,
                              if (!is.null(.endpt)) .endpt$sclIdx))
  .fixIdx <- which(.thetas$fix) - 1L
  # mixture component typical values (tka_k) and probabilities (p_k) are estimated
  # by the mixture M-step, not the structural (beta) re-solve.
  .mixIdx <- if (is.null(.mix)) integer(0) else c(.mix$muCompIdx, .mix$probIdx)
  structIdx <- setdiff(seq_len(nTheta) - 1L,
                       c(muIdx, covCoefIdx, as.integer(.resIdx), as.integer(.fixIdx),
                         as.integer(.mixIdx)))
  structIdx <- as.integer(structIdx)
  list(base = base, nTheta = nTheta, nEta = nEta, etaIdx = etaIdx, omega0 = omega0,
       muIdx = muIdx, mu0 = .thetas$est[muIdx + 1L],
       addSdIdx = addSdIdx, addSd0 = addSd0, errType = errType,
       propSdIdx = propSdIdx, propSd0 = propSd0, errName = errName,
       lambdaIdx = lambdaIdx, lambda0 = lambda0,
       tbsYj = .tbs$yj, tbsLow = .tbs$low, tbsHi = .tbs$hi,
       powIdx = powIdx, pow0 = pow0, endpt = .endpt,
       structIdx = structIdx, struct0 = .thetas$est[structIdx + 1L],
       fixIdx = as.integer(.fixIdx), fixNames = .thetas$name[.thetas$fix],
       covCoefNames = covCoefNames, covNames = covNames, covCoefIdx = covCoefIdx,
       covCoef0 = if (length(covCoefIdx)) .thetas$est[covCoefIdx + 1L] else numeric(0),
       thetaNames = .thetas$name, etaNames = .etas$name, muNames = .muName,
       mix = .mix)
}

#' Fit a mu-referenced model with RPEM (K=1 core).
#'
#' @param ui rxode2 UI object.
#' @param data Data frame (with a DV column).
#' @param control `rpemControl()`.
#' @return list of estimates (`mu`, `omega`, `addSd`) plus per-iteration traces.
#' @noRd
#' Form the RPEM Fisher-score covariance from the per-subject score matrix.
#'
#' The empirical Fisher information is `I = t(S) %*% S` and the covariance its
#' inverse.  `colNames` names the columns of `S` in order (typical value / covariate
#' coefs + residual sd + `om.<eta>` variances); the theta-named rows drive the
#' parFixedDf SEs and the `om.<eta>` rows carry the omega-variance SEs (reported in
#' `$cov`, as SAEM does).  NULL if the information is not invertible / the covariance
#' is not a valid (finite, non-negative-diagonal) matrix.
#' @noRd
.rpemFisherCov <- function(S, colNames) {
  .I <- crossprod(S)
  .cov <- tryCatch(solve(.I), error = function(e) NULL)
  if (is.null(.cov) || anyNA(.cov) || any(!is.finite(.cov)) || any(diag(.cov) < 0))
    return(NULL)
  dimnames(.cov) <- list(colNames, colNames)
  list(cov = .cov)
}

.rpemFit <- function(ui, data, control = rpemControl()) {
  .cl <- .rpemClassify(ui)
  if (!is.null(.cl$mix)) return(.rpemFitMix(ui, data, .cl, control))
  .m <- ui$rpemRxModel$predOnly
  .nm <- c(paste0("THETA[", seq_len(.cl$nTheta), "]"),
           paste0("ETA[", seq_len(.cl$nEta), "]"))
  .e <- new.env()
  .e$predOnly <- .m
  .e$rxControl <- rxode2::rxControl(atol = control$atol, rtol = control$rtol,
                                    cores = control$cores)
  .e$param <- stats::setNames(.cl$base, .nm)
  .e$data <- data
  .idColF <- names(data)[tolower(names(data)) == "id"][1]
  if (is.na(.idColF)) stop("data must have an ID column")
  .nsub <- length(unique(data[[.idColF]]))
  # censoring (BLQ): any non-zero CENS column -> the E-step uses the censored
  # likelihood (doCensNormal1) and the residual M-step maximizes it.
  .censColF <- names(data)[tolower(names(data)) == "cens"][1]
  .hasCens <- !is.na(.censColF) && any(data[[.censColF]] != 0, na.rm = TRUE)
  # Draw the E-step etas in R (deterministic given the seed) so the sampling RNG is
  # independent of the solve's core count -- the fit is reproducible across cores.
  # Optional mode-centered importance sampling (impInflate > 0): draw z ~ N(0,
  # impInflate*Omega) then shift each subject's block by its EBE (posterior mean from
  # the previous iteration).  impInflate == 0 keeps the paper's prior sampling
  # (.cInf = 1, ebe = 0 -> draw is exactly N(0, Omega)).
  .modeIS <- control$impInflate > 0
  .cInf <- if (.modeIS) control$impInflate else 1
  .drawEtas <- function(omega, ebe) {
    .z <- matrix(rxode2::rxRmvn(.nsub * control$nGauss, mu = rep(0, .cl$nEta),
                                sigma = .cInf * omega), ncol = .cl$nEta)
    .z + ebe[rep(seq_len(.nsub), each = control$nGauss), , drop = FALSE]
  }

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
  if (.multi) {
    .endptIdx <- .rpemEndptIndex(data, .cl); sdVec <- .cl$endpt$scl0
    propVec <- .cl$endpt$prop0
    .combE <- .cl$endpt$errType %in% c(2L, 3L, 4L)  # endpoints with a 2nd residual param
  }
  .structOn <- length(.cl$structIdx) > 0L
  niter <- control$niter
  coefTr <- if (.useReg) matrix(0, niter, ncol(.design)) else NULL
  muTr <- matrix(0, niter, .cl$nEta); omTr <- matrix(0, niter, .cl$nEta)
  sdTr <- numeric(niter); propTr <- numeric(niter); lamTr <- numeric(niter)
  powTr <- numeric(niter); llTr <- numeric(niter)
  sdMat <- if (.multi) matrix(0, niter, .cl$endpt$nEndpt) else NULL
  propMat <- if (.multi) matrix(NA_real_, niter, .cl$endpt$nEndpt) else NULL
  betaMat <- if (.structOn) matrix(0, niter, length(.cl$structIdx)) else NULL
  ebe <- matrix(0, .nsub, .cl$nEta)     # IS proposal center; updated each E-step when on
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
    if (.multi) {
      base[.cl$endpt$sclIdx + 1L] <- sdVec
      if (any(.combE)) base[.cl$endpt$propIdx[.combE] + 1L] <- propVec[.combE]
    }
    rxode2::rxSetSeed(control$seed + .it)
    .etaMat <- .drawEtas(omega, ebe)
    .est <- rpemEstepK1Draw(.e, base, .cl$etaIdx, .etaMat, control$nGauss, control$cores,
                            ebe, diag(as.matrix(omega)), .cInf)
    if (.modeIS) ebe <- .est$ebe        # posterior mean -> next proposal center
    # numeric M-step for non-mu-ref structural fixed effects, while the E-step
    # solve is still loaded (before the MH step's rxRmvn draw clobbers it).
    if (.structOn) {
      .bt <- rpemMstepBeta(base, .cl$etaIdx, .cl$structIdx, base[.cl$structIdx + 1L])
      base[.cl$structIdx + 1L] <- .bt; betaMat[.it, ] <- .bt
    }
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
      .prop0 <- ifelse(is.na(propVec), 0.0, propVec)
      .ms <- rpemMstepK1Multi(.design, coefs, .endptIdx, .cl$endpt$errType,
                              sdVec, .prop0, control$nMH, control$mhBurn)
      coefs <- .ms$coefs; omega <- matrix(.ms$omega, 1, 1)
      sdVec <- .ms$sd; propVec <- .ms$propSd
      coefTr[.it, ] <- coefs; muTr[.it, ] <- coefs[1]; omTr[.it, ] <- .ms$omega
      sdMat[.it, ] <- sdVec; propMat[.it, ] <- propVec
    } else if (.useReg) {
      .ms <- if (.hasCens && .cl$errType %in% c(0L, 1L))
        rpemMstepK1Cens(.design, coefs, .cl$errType, addSd, control$nMH, control$mhBurn)
      else rpemMstepK1Reg(.design, coefs, .cl$errType, control$nMH, control$mhBurn)
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
  endptPropHat <- if (.multi)
    apply(propMat[.w, , drop = FALSE], 2, function(x) if (all(is.na(x))) NA_real_ else mean(x)) else NULL
  structHat <- if (.structOn) colMeans(betaMat[.w, , drop = FALSE]) else NULL
  covCoefHat <- if (.useReg && length(.cl$covCoefIdx))
                  colMeans(coefTr[.w, -1, drop = FALSE]) else numeric(0)

  # One final E-step at the converged estimates to compute per-subject EBEs
  # (posterior-mean etas, Eq 53): EBE_i = sum_j eta_ij * w_ij with the
  # self-normalized importance weights w_ij = softmax_j(log p_ij).
  base[.cl$muIdx + 1L] <- muHat
  if (.structOn) base[.cl$structIdx + 1L] <- structHat
  if (length(covCoefHat)) base[.cl$covCoefIdx + 1L] <- covCoefHat
  if (!.multi) base[.cl$addSdIdx + 1L] <- sdHat
  if (.comb) base[.cl$propSdIdx + 1L] <- propHat
  if (.tbs) base[.cl$lambdaIdx + 1L] <- lambdaHat
  if (.pow) base[.cl$powIdx + 1L] <- powerHat
  if (.multi) {
    base[.cl$endpt$sclIdx + 1L] <- endptSdHat
    if (any(.combE)) base[.cl$endpt$propIdx[.combE] + 1L] <- endptPropHat[.combE]
  }
  omegaHat <- diag(omHat, .cl$nEta)
  rxode2::rxSetSeed(control$seed)
  .feEta <- .drawEtas(omegaHat, ebe)
  .fe <- rpemEstepK1Draw(.e, base, .cl$etaIdx, .feEta, control$nGauss, control$cores,
                         ebe, diag(as.matrix(omegaHat)), .cInf)
  # Fisher-score covariance (design/rpem/08) from the converged samples, computed
  # while they are still loaded (before rpemFree).  Analytic (importance-weighted
  # complete-data) scores for mu / diagonal Omega / residual (add, prop, combined,
  # power) via rpemFisherReg (single eta + covariates) or rpemFisherDiag (multi eta).
  # Non-mu-ref structural betas have no analytic score here (no stored dcp/dbeta), so
  # their per-subject marginal-loglik score is appended by a common-random-number
  # finite difference (reuse .feEta; only the beta moves).  TBS residual and mixtures
  # still keep the FOCEI-covariance SEs.
  .fisher <- NULL
  if (.cl$errType %in% c(0L, 1L, 2L, 4L)) {
    .omNames <- paste0("om.", .cl$etaNames)
    # residual parameter(s) + their theta names, ordered to match the C++ score
    # columns: combined = (add.sd, prop.sd); power = (prop.sd, power); else single sd
    if (.cl$errType == 2L) {
      .resPar <- c(sdHat, propHat)
      .resNm <- .cl$thetaNames[c(.cl$addSdIdx, .cl$propSdIdx) + 1L]
    } else if (.cl$errType == 4L) {
      .resPar <- c(sdHat, powerHat)
      .resNm <- .cl$thetaNames[c(.cl$addSdIdx, .cl$powIdx) + 1L]
    } else {
      .resPar <- sdHat
      .resNm <- .cl$thetaNames[.cl$addSdIdx + 1L]
    }
    if (.useReg) {                                    # nEta == 1 (+ optional covariates)
      .S <- rpemFisherReg(.design, c(muHat, covCoefHat), omHat[1], .cl$errType, .resPar)
      .colN <- c(.cl$muNames, .cl$covCoefNames, .resNm, .omNames)
    } else {                                          # nEta > 1, diagonal Omega
      .S <- rpemFisherDiag(muHat, omHat, .cl$errType, .resPar)
      .colN <- c(.cl$muNames, .resNm, .omNames)
    }
    if (.structOn) {
      # per-subject marginal loglik s_i = logsumexp_j logp_ij (the -log(nG) constant
      # cancels in the difference); score_i,m = (s_i(beta_m + h) - s_i(beta_m))/h with
      # the same drawn etas .feEta, so only beta_m changes between the two solves.
      .nGf <- control$nGauss; .nsubF <- length(.fe$logp) / .nGf
      .lseSub <- function(lp) vapply(seq_len(.nsubF), function(.i) {
        .v <- lp[((.i - 1L) * .nGf + 1L):(.i * .nGf)]; .mx <- max(.v)
        .mx + log(sum(exp(.v - .mx)))
      }, numeric(1))
      .s0 <- .lseSub(.fe$logp)
      .Sb <- matrix(0, .nsubF, length(.cl$structIdx))
      for (.m in seq_along(.cl$structIdx)) {
        .h <- 1e-4 * max(abs(structHat[.m]), 1)
        .bp <- base; .bp[.cl$structIdx[.m] + 1L] <- .bp[.cl$structIdx[.m] + 1L] + .h
        .ep <- rpemEstepK1Draw(.e, .bp, .cl$etaIdx, .feEta, control$nGauss, control$cores,
                               ebe, diag(as.matrix(omegaHat)), .cInf)
        .Sb[, .m] <- (.lseSub(.ep$logp) - .s0) / .h
      }
      .S <- cbind(.S, .Sb)
      .colN <- c(.colN, .cl$thetaNames[.cl$structIdx + 1L])
    }
    .fisher <- .rpemFisherCov(.S, .colN)
  }
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
       endptSd = endptSdHat, endptProp = endptPropHat,
       struct = if (.structOn) stats::setNames(structHat, .cl$thetaNames[.cl$structIdx + 1L]) else NULL,
       covCoef = stats::setNames(covCoefHat, .cl$covCoefNames),
       ebe = ebe, fisher = .fisher,
       lnL = llTr, muTrace = muTr, omegaTrace = omTr, sdTrace = sdTr,
       classify = .cl)
}

#' Fit a mix() split-ETA mixture with RPEM.
#'
#' Mixture EM (paper section on finite mixtures): the E-step solves each Monte
#' Carlo sample once per component (rpemEstepMixDraw), and the M-step (rpemMstepMix)
#' runs a Metropolis-Hastings over the joint (subject, sample, component) posterior
#' to update the per-component typical values (mu_k), the shared BSV (Omega), the
#' mixture weights (w_k) and the shared residual.  Requires a single mixed random
#' effect, an additive / proportional / lognormal residual, and no non-mixture
#' structural fixed effects to estimate (fix() any structural typical values).
#' @noRd
.rpemFitMix <- function(ui, data, .cl, control = rpemControl()) {
  .mix <- .cl$mix; K <- .mix$K
  if (!.cl$errType %in% c(0L, 1L, 6L))
    stop("RPEM mixtures currently support additive, proportional or lognormal residuals")
  if (length(.cl$structIdx) > 0L)
    stop("RPEM mixtures require non-mixture structural typical values to be fix()ed")
  .m <- ui$rpemRxModel$predOnly
  .nm <- c(paste0("THETA[", seq_len(.cl$nTheta), "]"),
           paste0("ETA[", seq_len(.cl$nEta), "]"))
  .e <- new.env()
  .e$predOnly <- .m
  .e$rxControl <- rxode2::rxControl(atol = control$atol, rtol = control$rtol,
                                    cores = control$cores)
  .e$param <- stats::setNames(.cl$base, .nm)
  .e$data <- data
  .idColF <- names(data)[tolower(names(data)) == "id"][1]
  if (is.na(.idColF)) stop("data must have an ID column")
  .nsub <- length(unique(data[[.idColF]]))
  .drawEtas <- function(omega) rxode2::rxRmvn(.nsub * control$nGauss, mu = 0, sigma = omega)

  base <- .cl$base; muK <- .mix$muComp0; w <- .mix$w0
  omega <- matrix(.cl$omega0[1, 1], 1, 1); addSd <- .cl$addSd0
  niter <- control$niter
  muMat <- matrix(0, niter, K); wMat <- matrix(0, niter, K)
  omTr <- numeric(niter); sdTr <- numeric(niter); llTr <- numeric(niter)
  for (.it in seq_len(niter)) {
    base[.mix$muCompIdx + 1L] <- muK
    base[.cl$addSdIdx + 1L] <- addSd
    rxode2::rxSetSeed(control$seed + .it)
    .etaMat <- .drawEtas(omega)
    .est <- rpemEstepMixDraw(.e, base, .cl$etaIdx, .etaMat, control$nGauss, control$cores, K, w)
    .ms <- rpemMstepMix(muK, w, .cl$errType, control$nMH, control$mhBurn)
    muK <- .ms$muK; omega <- matrix(.ms$omega, 1, 1); w <- .ms$w; addSd <- .ms$addSd
    muMat[.it, ] <- muK; wMat[.it, ] <- w; omTr[.it] <- .ms$omega
    sdTr[.it] <- addSd; llTr[.it] <- .est$lnL
  }
  .k <- min(control$collect, niter); .wi <- (niter - .k + 1L):niter
  muHat <- colMeans(muMat[.wi, , drop = FALSE])
  wHat <- colMeans(wMat[.wi, , drop = FALSE])
  omHat <- mean(omTr[.wi]); sdHat <- mean(sdTr[.wi])

  # Final E-step at the converged estimates: EBEs + per-subject component posteriors.
  base[.mix$muCompIdx + 1L] <- muHat; base[.cl$addSdIdx + 1L] <- sdHat
  omegaHat <- matrix(omHat, 1, 1)
  rxode2::rxSetSeed(control$seed)
  .feEta <- .drawEtas(omegaHat)
  .fe <- rpemEstepMixDraw(.e, base, .cl$etaIdx, .feEta, control$nGauss, control$cores, K, wHat)
  rpemFree()
  # component posteriors tau_ik = w_k n_ik / sum_k, mixNum = argmax.
  .lw <- sweep(.fe$lognik, 2, log(wHat), "+")
  .tau <- exp(.lw - apply(.lw, 1, max)); .tau <- .tau / rowSums(.tau)
  mixNum <- max.col(.tau, ties.method = "first")
  # Per-component EBEs (self-normalized importance weights within each component):
  # EBE_ik = sum_j eta_ij * softmax_j(log p(Y_i | k, eta_ij)).  ebeK is nsub x K.
  .nG <- control$nGauss
  .etaVec <- as.numeric(.fe$eta)
  ebeK <- matrix(0, .nsub, K)
  for (.i in seq_len(.nsub)) {
    .idx <- ((.i - 1L) * .nG + 1L):(.i * .nG)
    for (.kc in seq_len(K)) {
      .lp <- .fe$logp[.idx, .kc]
      .wt <- exp(.lp - max(.lp)); .wt <- .wt / sum(.wt)
      ebeK[.i, .kc] <- sum(.etaVec[.idx] * .wt)
    }
  }
  # Display EBE = winning component's; FOCEI eval etaMat = per-component EBEs stacked
  # (component 1 for all subjects, then component 2, ...): (K*nsub) x neta.
  ebe <- matrix(ebeK[cbind(seq_len(.nsub), mixNum)], .nsub, 1L,
                dimnames = list(NULL, .cl$etaNames))
  ebeStack <- matrix(as.numeric(ebeK), K * .nsub, 1L, dimnames = list(NULL, .cl$etaNames))

  list(mu = stats::setNames(muHat[1], .cl$muNames),   # nominal (comp-1) mu for the eta
       omega = stats::setNames(omHat, .cl$etaNames),
       addSd = sdHat, propSd = if (.cl$errType == 1L) sdHat else NA_real_,
       lambda = NA_real_, power = NA_real_,
       endptSd = NULL, endptProp = NULL, struct = NULL,
       covCoef = stats::setNames(numeric(0), character(0)), ebe = ebe,
       lnL = llTr, muTrace = muMat, omegaTrace = matrix(omTr, ncol = 1), sdTrace = sdTr,
       classify = .cl,
       mix = list(K = K, muK = stats::setNames(muHat, .mix$muNames),
                  w = wHat, probNames = .mix$probNames, mixNum = mixNum, tau = .tau,
                  ebeK = ebeK, ebeStack = ebeStack))
}

#' Populate the mixture fit fields (mixList / mixNum / probabilities / iCov) on
#' the RPEM fit env from the RPEM posterior component weights (tau) and mixture
#' probabilities (w).  Analogous to `.saemMixFix` for the single-mixed-eta case
#' (no eta-group collapsing needed since the components share one random effect).
#' @noRd
.rpemMixSetFit <- function(env, .cl, rfit) {
  .mix <- rfit$mix; .K <- .mix$K
  .etaNames <- .cl$etaNames
  .tau <- .mix$tau                         # nsub x K posterior weights
  .best <- .mix$mixNum
  # subject IDs + EBEs from the fit ranef (post-eval) or the pre-eval etaObf.
  .src <- if (!is.null(env$ranef)) as.data.frame(env$ranef) else as.data.frame(env$etaObf)
  .ids <- .src$ID
  .etaDf <- .src[, .etaNames, drop = FALSE]
  env$mixProbabilities <- unname(.mix$w)
  # ranef gains a mixnum column; mixList: one frame per component (ID, etas, prob).
  .ranef <- cbind(data.frame(ID = .ids), .etaDf); .ranef$mixnum <- as.integer(.best)
  assign("ranef", .ranef, envir = env)
  .mixList <- lapply(seq_len(.K), function(k) {
    .r <- cbind(data.frame(ID = .ids), .etaDf, data.frame(prob = .tau[, k]))
    names(.r) <- c("ID", .etaNames, "prob"); row.names(.r) <- NULL; .r
  })
  names(.mixList) <- paste0("mix", seq_len(.K))
  assign("mixList", .mixList, envir = env)
  .mixNum <- data.frame(ID = .ids, mixnum = as.integer(.best)); row.names(.mixNum) <- NULL
  assign("mixNum", .mixNum, envir = env)
  # expected etas (for shrinkage): sum_k prob_k * eta_k -- one shared eta, so the
  # posterior-weighted eta reduces to the EBE already in etaObf.
  .etaExp <- cbind(data.frame(ID = .ids), .etaDf)
  assign("etaExpected", .etaExp, envir = env)
  # iCov fixes each subject's mixture component during the eval table solve.
  assign("mixIcov", data.frame(ID = as.integer(.ids), mixest = as.integer(.best)),
         envir = env)
  invisible(NULL)
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
  # Mixtures produce IPRED/PRED tables but skip CWRES/NPDE (as SAEM does): those add
  # a FOCEi re-fit whose per-subject mixture etaMat carries a mixnum column that the
  # inner solver rejects.  Disable them in the table control for a mixture fit.
  if (!is.null(.cl$mix) && inherits(.ret$table, "tableControl")) {
    .ret$table$cwres <- FALSE
    .ret$table$npde <- FALSE
  }
  # full theta (named over all theta names): mu-referenced -> RPEM mu, additive
  # residual -> RPEM add.sd, held structural -> ini values.
  .tn <- .cl$thetaNames
  .ft <- stats::setNames(.cl$base[seq_along(.tn)], .tn)
  .ft[.cl$muNames] <- rfit$mu
  if (.cl$errType == 5L) {
    .ft[.tn[.cl$endpt$sclIdx + 1L]] <- rfit$endptSd
    .cE <- which(.cl$endpt$errType %in% c(2L, 3L, 4L))  # 2nd param: prop.sd (combined) or exponent (power)
    if (length(.cE)) .ft[.tn[.cl$endpt$propIdx[.cE] + 1L]] <- rfit$endptProp[.cE]
  } else .ft[.tn[.cl$addSdIdx + 1L]] <- rfit$addSd
  if (.cl$errType == 2L) .ft[.tn[.cl$propSdIdx + 1L]] <- rfit$propSd
  if (.cl$errType == 3L) .ft[.tn[.cl$lambdaIdx + 1L]] <- rfit$lambda
  # power: exponent set here; the scale (prop.sd) is set via the addSd slot above.
  if (.cl$errType == 4L) .ft[.tn[.cl$powIdx + 1L]] <- rfit$power
  if (length(rfit$struct)) .ft[.tn[.cl$structIdx + 1L]] <- rfit$struct
  if (length(rfit$covCoef)) .ft[.cl$covCoefNames] <- rfit$covCoef
  # mixture: per-component typical values (tka_k) and mixture probabilities (p_k).
  if (!is.null(.cl$mix)) {
    .ft[.cl$mix$muNames] <- rfit$mix$muK
    # the eval reads fullTheta[p1] through its mlogit parameter transform before the
    # solve (yet validates it as natural at setup); mexpit(w) satisfies both, since
    # mlogit(mexpit(w)) == w recovers the natural weight for the model.
    if (length(.cl$mix$probNames))
      .ft[.cl$mix$probNames] <- rxode2::mexpit(rfit$mix$w[seq_along(.cl$mix$probNames)])
  }
  # Seed the ui iniDf with the RPEM estimates so the eval-only FOCEI (maxOuter=0)
  # reports them (it starts from -- and holds -- the iniDf theta).  Only user-fixed
  # (fix()) thetas are held; mark them fixed so they report as held (NA SE).
  .heldNames <- .cl$fixNames
  .uiD <- rxode2::rxUiDecompress(ui)
  .idf <- .uiD$iniDf
  .thRow <- !is.na(.idf$ntheta)
  .idf$est[.thRow] <- .ft[.idf$name[.thRow]]
  if (length(.heldNames) > 0L) .idf$fix[.idf$name %in% .heldNames] <- TRUE
  assign("iniDf", .idf, envir = .uiD)
  ui <- rxode2::rxUiCompress(.uiD)
  .ret$ui <- ui
  .ret$fullTheta <- .ft
  # omega (diagonal) named over etas
  .om <- diag(rfit$omega, .cl$nEta)
  dimnames(.om) <- list(.cl$etaNames, .cl$etaNames)
  .ret$omega <- .om
  # per-subject EBEs -> etaMat for the FOCEI eval (nsub x neta; display EBE = winning
  # component for a mixture).  FOCEI expands this to K*nsub for the mixture itself.
  .isMix <- !is.null(.cl$mix)
  .eb <- rfit$ebe
  colnames(.eb) <- .cl$etaNames
  .ret$.etaMat <- .eb
  .ret$.etaMatBase <- .eb
  .ret$etaObf <- data.frame(ID = seq_len(nrow(.eb)),
                            stats::setNames(as.data.frame(.eb), .cl$etaNames),
                            OBJI = NA)
  # mixture: build mixList / mixNum / mixIcov / probabilities BEFORE the eval (mirrors
  # SAEM's .saemMixFix), then hand the FOCEI eval an N-row etaMat via the control only
  # (remove env$.etaMat so no mixnum column leaks into it, exactly as SAEM does).
  if (.isMix) .rpemMixSetFit(.ret, .cl, rfit)
  # RPEM Fisher-score covariance (design/rpem/08): when available, install it as the
  # fit's cov (the paper's native SEs) and skip the FOCEI covariance calc; otherwise
  # covMethod "r,s" gives FOCEI-covariance SEs at the fixed RPEM estimates.  Mixtures
  # skip the covariance entirely (perturbing the mixture probability breaks the solve),
  # matching SAEM's covMethod=0 for the eval.
  # The FOCEI eval runs with covMethod "r,s" for non-mixtures (it also settles the
  # displayed estimates); when Fisher-score SEs are available they OVERWRITE the
  # FOCEI covariance / parFixedDf SEs afterward (.rpemInstallFisherCov).
  .fisherOn <- !is.null(rfit$fisher)
  # eventSens="jump" enables rxode2's analytic dose-parameter ("jump") sensitivities
  # for the FOCEI eval's gradient model (residual tables + R/S covariance) -- set
  # explicitly so dose-parameter models (bioavailability / lag / rate / dur) get the
  # correct event sensitivities regardless of the foceiControl default.
  .foceiControl <- foceiControl(maxOuterIterations = 0L, maxInnerIterations = 0L,
                                covMethod = if (.isMix) 0L else "r,s",
                                etaMat = .eb, scaleTo = 0,
                                skipCov = if (.isMix) ui$foceiSkipCov else NULL,
                                calcTables = TRUE, interaction = 1L, eventSens = "jump",
                                rxControl = .rxControl, est = "rpem")
  if (.isMix) {
    # the FOCEI eval must read the etaMat from the control, not env$.etaMat.
    if (exists(".etaMat", envir = .ret, inherits = FALSE)) rm(".etaMat", envir = .ret)
    if (exists(".etaMatBase", envir = .ret, inherits = FALSE)) rm(".etaMatBase", envir = .ret)
  }
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
  # Install the RPEM Fisher-score covariance AFTER the fit is built (so the eval-only
  # finalize does not clobber the estimates) and patch the parFixedDf SE/%RSE/CI from
  # it -- mirrors SAEM's .saemInstallFullCov.
  if (.fisherOn) .rpemInstallFisherCov(.out, rfit$fisher$cov)
  .out
}

#' Install the RPEM Fisher-score covariance on a built fit and surface its SEs in
#' parFixedDf (SE / %RSE / CI) for the theta rows it covers.  Mirrors
#' `.saemInstallFullCov`; `cov` is the theta-block covariance (typical value +
#' covariate coefs + residual sd) named over the estimated thetas.
#' @noRd
.rpemInstallFisherCov <- function(fit, cov) {
  .env <- if (rxode2::rxIs(fit, "nlmixr2FitData")) fit$env else fit
  if (!is.environment(.env) || is.null(cov)) return(invisible())
  .env$cov <- cov
  .env$covMethod <- "fisher"
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
