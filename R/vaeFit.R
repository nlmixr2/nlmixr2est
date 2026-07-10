# vaeFit.R -- VAE training: encoder init, one ELBO gradient step (encoder C++
# forward/backward + rxode2 decoder), Adam, closed-form M-step, and the
# burn-in -> KL-anneal -> smoothing schedule. Orchestration is thin R; the heavy
# work is the C++ encoder BPTT (vaeEncoderFwdBwd) and the rxode2 decoder solve.

#' Initialize encoder parameters (RNG seeded once by the caller under rxWithSeed)
#' @noRd
.vaeEncoderInitParams <- function(zDim, hDim, nCov, zPop, sigma0) {
  nOff <- as.integer(zDim * (zDim - 1L) / 2L)
  outDim <- 2L * zDim + nOff
  .sd <- 1 / sqrt(hDim)
  list(
    Wih = matrix(stats::rnorm(4L * hDim * 2L, 0, .sd), 4L * hDim, 2L),
    Whh = matrix(stats::rnorm(4L * hDim * hDim, 0, .sd), 4L * hDim, hDim),
    bih = numeric(4L * hDim),
    bhh = numeric(4L * hDim),
    fcW = matrix(stats::rnorm(outDim * (hDim + nCov), 0, 1e-2), outDim, hDim + nCov),
    fcB = c(zPop, log(sigma0), numeric(nOff))
  )
}

#' Clamp a parameter vector to [lower, upper] (elementwise; NULL bounds = no-op).
#' @noRd
.vaeClamp <- function(v, lower, upper) {
  if (!is.null(lower)) v <- pmax(v, lower)
  if (!is.null(upper)) v <- pmin(v, upper)
  v
}

#' Assemble the full theta vector from current z_pop (structural) + a (residual)
#' @noRd
.vaeBuildTh <- function(prep, zPop, a) {
  th <- prep$th
  ## mixture etas (NA index) keep their fixed component thetas in `th`
  .ok <- !is.na(prep$zPopThetaIdx)
  th[prep$zPopThetaIdx[.ok]] <- zPop[.ok]
  if (length(prep$errThetaIdx) > 0L) th[prep$errThetaIdx] <- a
  th
}

#' One ELBO evaluation: encoder forward, decoder solve, KL, encoder backward.
#' Returns the loss, encoder param gradients, and quantities the M-step needs.
#' @param params encoder parameter list
#' @param prep .vaeDataPrep output
#' @param am decoder augmented model (.vaeDecoderModel)
#' @param zPop,omega,a current population parameters
#' @param alphaKL KL weight
#' @param eps fixed reparam noise [N, zDim]
#' @param withGrad compute the encoder backward (default TRUE)
#' @noRd
.vaeElboStep <- function(params, prep, am, zPop, omega, a, alphaKL, eps, withGrad = TRUE) {
  N <- prep$N; zDim <- prep$zDim
  ## zPop may be a shared vector (no covariates) or an N x zDim matrix of
  ## subject-specific KL centers (covariate model). The decoder baseline
  ## (THETA) is immaterial to f (only THETA+ETA=z matters), so use the intercept.
  if (is.matrix(zPop)) { zPopMat <- zPop; baseline <- colMeans(zPop) }
  else { zPopMat <- matrix(zPop, N, zDim, byrow = TRUE); baseline <- zPop }
  th <- .vaeBuildTh(prep, baseline, a)
  z0 <- matrix(0, N, zDim)
  ## forward (grads ignored)
  fw <- vaeEncoderFwdBwd(prep$dataIn, prep$lengths, prep$covIn, eps,
                         params$Wih, params$Whh, params$bih, params$bhh,
                         params$fcW, params$fcB, zDim, z0, z0)
  mu <- fw$mu; logSigma <- fw$logSigma; L <- fw$L; z <- fw$z

  ## decoder per subject -> p(x|z) and d(p_x_z)/dz
  gZdec <- matrix(0, N, zDim)
  pxz <- 0
  preds <- vector("list", N)
  for (i in seq_len(N)) {
    s <- prep$subj[[i]]
    eta <- z[i, ] - baseline
    E <- .vaeDecoderSolveSubject(am, th, eta, s$ev, s$times)
    if (is.null(E)) return(NULL)
    px <- .vaeDecoderPxz(E, s$y)
    pxz <- pxz + px$pxz
    gZdec[i, ] <- px$gEta
    preds[[i]] <- E$f
  }

  ## KL: p_z - q_z, centered on the (possibly subject-specific) z_pop
  ln2pi <- log(2 * pi)
  pz <- 0; qz <- 0
  for (i in seq_len(N)) {
    pz <- pz + 0.5 * sum((z[i, ] - zPopMat[i, ])^2 / omega + log(omega) + ln2pi)
    qz <- qz + 0.5 * sum(eps[i, ]^2 + ln2pi + 2 * logSigma[i, ])
  }
  DKL <- pz - qz
  loss <- pxz + alphaKL * DKL

  grads <- NULL
  if (withGrad) {
    gZ <- gZdec
    for (i in seq_len(N)) gZ[i, ] <- gZ[i, ] + alphaKL * (z[i, ] - zPopMat[i, ]) / omega
    gLS <- matrix(-alphaKL, N, zDim)
    bw <- vaeEncoderFwdBwd(prep$dataIn, prep$lengths, prep$covIn, eps,
                           params$Wih, params$Whh, params$bih, params$bhh,
                           params$fcW, params$fcB, zDim, gZ, gLS)
    grads <- list(Wih = bw$gWih, Whh = bw$gWhh, bih = bw$gbih, bhh = bw$gbhh,
                  fcW = bw$gFcW, fcB = bw$gFcB)
  }
  list(loss = loss, pxz = pxz, DKL = DKL, grads = grads,
       mu = mu, L = L, z = z, preds = preds)
}

#' Adam optimizer state (per encoder parameter block)
#' @noRd
.vaeAdamInit <- function(params) {
  lapply(params, function(p) {
    z <- p; z[] <- 0
    list(m = z, v = z)
  })
}

#' One Adam update over all encoder parameter blocks
#' @noRd
.vaeAdamStep <- function(params, grads, state, lr, t, b1 = 0.9, b2 = 0.999, eps = 1e-8) {
  for (nm in names(params)) {
    g <- grads[[nm]]
    state[[nm]]$m <- b1 * state[[nm]]$m + (1 - b1) * g
    state[[nm]]$v <- b2 * state[[nm]]$v + (1 - b2) * g^2
    mhat <- state[[nm]]$m / (1 - b1^t)
    vhat <- state[[nm]]$v / (1 - b2^t)
    params[[nm]] <- params[[nm]] - lr * mhat / (sqrt(vhat) + eps)
  }
  list(params = params, state = state)
}

#' Closed-form M-step (no covariate selection): z_pop = mean posterior mean;
#' omega_k = mean_i[(mu_ik - z_pop_k)^2 + (L_i L_i^T)_kk]; a = sqrt(SSR / Nobs).
#' `gamma` is the EMA smoothing rate (1 = direct update).
#' @noRd
.vaeMStep <- function(mu, L, preds, prep, zPop, omega, a, gamma) {
  N <- prep$N; zDim <- prep$zDim
  zPopCur <- colMeans(mu)
  omegaCur <- numeric(zDim)
  for (k in seq_len(zDim)) {
    v <- 0
    for (i in seq_len(N)) v <- v + (mu[i, k] - zPopCur[k])^2 + sum(L[k, , i]^2)
    omegaCur[k] <- v / N
  }
  aCur <- .vaeUpdateErr(preds, prep, a)
  ## non-mu-referenced (free) etas center at 0 (theta forced to 0); fixed omega
  ## entries are held at their current value (not estimated)
  if (!is.null(prep$isFree)) zPopCur[prep$isFree] <- 0
  if (!is.null(prep$omegaFix)) omegaCur[prep$omegaFix] <- omega[prep$omegaFix]
  if (!all(is.finite(zPopCur))) zPopCur <- zPop
  if (!all(is.finite(omegaCur))) omegaCur <- omega
  zPopCur <- .vaeClamp(zPopCur, prep$zPopLower, prep$zPopUpper)
  list(zPop = zPop + gamma * (zPopCur - zPop),
       omega = omega + gamma * (omegaCur - omega),
       a = a + gamma * (aCur - a))
}

#' Closed-form error-parameter M-step for additive / proportional / combined
#' residual models, robust to non-finite predictions (dropped, not poisoning the
#' estimate). Returns the error-param vector in `prep$errThetaIdx` order.
#'  add:      R = a^2               -> a = sqrt(mean(res^2))
#'  prop:     R = (b*f)^2           -> b = sqrt(mean((res/f)^2))
#'  combined: R = a^2 + (b*f)^2     -> nnls of res^2 on [1, f^2]
#' Other error types keep their current value (the inner likelihood still uses
#' them correctly; only their closed-form update is unavailable).
#' @noRd
.vaeUpdateErr <- function(preds, prep, a) {
  if (length(a) == 0L) return(a)
  res <- numeric(0); f <- numeric(0)
  for (i in seq_len(prep$N)) {
    r <- prep$subj[[i]]$y - preds[[i]]; ff <- preds[[i]]
    ok <- is.finite(r) & is.finite(ff)
    res <- c(res, r[ok]); f <- c(f, ff[ok])
  }
  if (length(res) == 0L) return(a)
  types <- prep$errType
  hasAdd <- which(types == "add"); hasProp <- which(types == "prop")
  aNew <- a
  if (length(hasAdd) && length(hasProp)) {
    ## combined: res^2 ~ a^2 + b^2 f^2 (non-negative least squares, 2 columns)
    X <- cbind(1, f^2); cf <- tryCatch(stats::lm.fit(X, res^2)$coefficients, error = function(e) c(NA, NA))
    v0 <- max(cf[1], .Machine$double.eps); v1 <- max(cf[2], .Machine$double.eps)
    if (is.finite(v0)) aNew[hasAdd[1]] <- sqrt(v0)
    if (is.finite(v1)) aNew[hasProp[1]] <- sqrt(v1)
  } else if (length(hasAdd)) {
    aNew[hasAdd[1]] <- sqrt(mean(res^2))
  } else if (length(hasProp)) {
    ok <- abs(f) > 1e-8
    if (any(ok)) aNew[hasProp[1]] <- sqrt(mean((res[ok] / f[ok])^2))
  }
  aNew[!is.finite(aNew)] <- a[!is.finite(aNew)]
  .vaeClamp(aNew, prep$errLower, prep$errUpper)
}

#' Closed-form M-step WITH BICc-ELBO covariate selection. For each parameter k,
#' enumerate the 2^nCov covariate subsets, fit OLS mu_ik ~ [1 | cov_S], and pick
#' the subset minimizing RSS_S/omega_k + log(N)*|S| (the per-parameter reduction
#' of the global BICc-ELBO L0 objective -- exact since the problem decouples by
#' parameter). Returns intercept + beta matrix + selected mask + subject-specific
#' z_pop matrix, and the residual-variance omega / a.
#' @noRd
.vaeMStepCov <- function(mu, L, preds, prep, omega, a, gamma) {
  N <- prep$N; zDim <- prep$zDim; nCov <- ncol(prep$covMat)
  X <- cbind(1, prep$covMat)                       # [N, 1 + nCov]
  intercept <- numeric(zDim); beta <- matrix(0, zDim, nCov)
  selected <- matrix(FALSE, zDim, nCov)
  zPopMat <- matrix(0, N, zDim)
  logN <- log(N)
  for (k in seq_len(zDim)) {
    ## a non-mu-referenced (free) eta is a theta=0-centered random effect: no
    ## intercept, no covariates (its structure is carried elsewhere in the model)
    if (!is.null(prep$isFree) && prep$isFree[k]) next
    yk <- mu[, k]; best <- NULL; bestScore <- Inf
    for (mask in 0:(2^nCov - 1L)) {
      S <- which(bitwAnd(mask, bitwShiftL(1L, seq_len(nCov) - 1L)) != 0L)
      cols <- c(1L, 1L + S)
      Xs <- X[, cols, drop = FALSE]
      fit <- stats::lm.fit(Xs, yk)
      score <- sum(fit$residuals^2) / omega[k] + logN * length(S)
      if (score < bestScore) { bestScore <- score; best <- list(S = S, coef = fit$coefficients, cols = cols) }
    }
    ## clamp the typical value (intercept) to its bounds
    intercept[k] <- .vaeClamp(best$coef[1], prep$zPopLower[k], prep$zPopUpper[k])
    best$coef[1] <- intercept[k]
    if (length(best$S) > 0L) { beta[k, best$S] <- best$coef[-1]; selected[k, best$S] <- TRUE }
    zPopMat[, k] <- X[, best$cols, drop = FALSE] %*% best$coef
  }
  omegaCur <- numeric(zDim)
  for (k in seq_len(zDim)) {
    v <- 0
    for (i in seq_len(N)) v <- v + (mu[i, k] - zPopMat[i, k])^2 + sum(L[k, , i]^2)
    omegaCur[k] <- v / N
  }
  aCur <- .vaeUpdateErr(preds, prep, a)
  if (!is.null(prep$omegaFix)) omegaCur[prep$omegaFix] <- omega[prep$omegaFix]
  if (!all(is.finite(omegaCur))) omegaCur <- omega
  list(intercept = intercept, beta = beta, selected = selected, zPopMat = zPopMat,
       omega = omega + gamma * (omegaCur - omega), a = a + gamma * (aCur - a))
}

#' Tracked population parameters (structural typical values, omega diagonal,
#' residual error) as a single named vector -- one parameter-history row.
#' @noRd
.vaeParRow <- function(zPop, omega, a, parInfo) {
  .z <- if (length(parInfo$structIdx)) setNames(zPop[parInfo$structIdx], parInfo$structNames) else numeric(0)
  c(.z, setNames(omega, parInfo$omegaNames), setNames(as.numeric(a), parInfo$aNames))
}

#' Train the VAE: burn-in (encoder-only, tiny KL) -> main EM (KL anneal + M-step)
#' -> EMA smoothing. Returns the trained encoder + population estimates + traces,
#' and the full parameter-history walk (`parHist`).
#' @noRd
.vaeTrain <- function(prep, innerEnv, control, nMix = 1L, mixProb = 1,
                      parInfo = NULL) {
  ## RNG is seeded ONCE for the whole estimation in nlmixr2Est.vae (rxWithSeed),
  ## which also covers the model's own random draws and restores the caller's seed
  zDim <- prep$zDim; hDim <- control$hiddenDim; nCov <- ncol(prep$covIn); N <- prep$N
  sigma0 <- if (is.null(control$sigma0)) rep(0.1, zDim) else rep_len(control$sigma0, zDim)
  params <- .vaeEncoderInitParams(zDim, hDim, nCov, prep$zPop, sigma0)
  adam <- .vaeAdamInit(params)
  zPop <- prep$zPop; omega <- prep$omega; a <- prep$a
  Lg <- control$nGradStep
  .t <- 0L; last <- NULL

  ## The parameter-history walk is ALWAYS captured (it is central to this method)
  ## via the shared iteration-print machinery (scale.h), so the walk prints like
  ## saem/focei and becomes standard parHistData. `parInfo` only supplies nicer
  ## structural names (the mu-referenced theta names) and the back-transform
  ## codes (`xform`, from .iterPrintXParFromUi) -- default to the eta names.
  if (is.null(parInfo)) {
    .sIdx <- which(!prep$isFree)
    parInfo <- list(structIdx = .sIdx, structNames = prep$etaNames[.sIdx],
                    omegaNames = paste0("o(", prep$etaNames, ")"), aNames = names(prep$a))
  }
  .row0 <- .vaeParRow(zPop, omega, a, parInfo)
  vaeIterPrintStart_(.row0, names(.row0), control$iterPrintControl, parInfo$xform)

  ## burn-in: encoder-only training with a tiny fixed KL weight
  for (it in seq_len(control$itersBurnIn)) {
    for (l in seq_len(Lg)) {
      eps <- matrix(stats::rnorm(N * zDim), N, zDim)
      st <- .vaeElboStepInner(params, prep, innerEnv, zPop, omega, a, 0.001, eps, control, nMix, mixProb)
      if (is.null(st)) stop("est=\"vae\" inner solve failed during burn-in", call. = FALSE)
      .t <- .t + 1L
      .ad <- .vaeAdamStep(params, st$grads, adam, control$burnInLearningRate, .t)
      params <- .ad$params; adam <- .ad$state
      .ms <- .vaeMStep(st$mu, st$L, st$preds, prep, zPop, omega, a, 1)
      zPop <- .ms$zPop; omega <- .ms$omega; a <- .ms$a
      last <- st
    }
    vaeIterPrintRow_(.vaeParRow(zPop, omega, a, parInfo), last$pxz + last$DKL, "Burn in")
  }

  ## main EM (optionally with BICc-ELBO covariate selection)
  doCov <- isTRUE(control$covariateSelection) && ncol(prep$covMat) > 0L
  intercept <- zPop; beta <- matrix(0, zDim, ncol(prep$covMat))
  selected <- matrix(FALSE, zDim, ncol(prep$covMat))
  zPopArg <- zPop
  nMain <- control$iters
  elboTrace <- numeric(nMain)
  for (it in seq_len(nMain)) {
    gamma <- if (it <= control$gammaIter) 1 else 1 / (1 + it - control$gammaIter)
    if (doCov) {
      .ms <- .vaeMStepCov(last$mu, last$L, last$preds, prep, omega, a, gamma)
      intercept <- .ms$intercept; beta <- .ms$beta; selected <- .ms$selected
      zPopArg <- .ms$zPopMat; zPop <- intercept; omega <- .ms$omega; a <- .ms$a
    } else {
      .ms <- .vaeMStep(last$mu, last$L, last$preds, prep, zPop, omega, a, gamma)
      zPop <- .ms$zPop; zPopArg <- zPop; omega <- .ms$omega; a <- .ms$a
    }
    alphaKL <- if (it <= control$klWarmup) 0.01 + 0.99 * (it - 1) / max(1, control$klWarmup - 1) else 1
    elbos <- numeric(Lg)
    for (l in seq_len(Lg)) {
      eps <- matrix(stats::rnorm(N * zDim), N, zDim)
      st <- .vaeElboStepInner(params, prep, innerEnv, zPopArg, omega, a, alphaKL, eps, control, nMix, mixProb)
      if (is.null(st)) stop("est=\"vae\" inner solve failed during training", call. = FALSE)
      .t <- .t + 1L
      .ad <- .vaeAdamStep(params, st$grads, adam, control$learningRate, .t)
      params <- .ad$params; adam <- .ad$state
      elbos[l] <- st$pxz + st$DKL
      last <- st
    }
    elboTrace[it] <- mean(elbos)
    .phase <- if (it <= control$klWarmup) {
      "KL anneal"
    } else if (it <= control$gammaIter) {
      "EM"
    } else {
      "Smooth"
    }
    vaeIterPrintRow_(.vaeParRow(zPop, omega, a, parInfo), elboTrace[it], .phase)
  }

  parHist <- vaeIterPrintGet_(isTRUE(control$print >= 1L))

  zPopMat <- if (is.matrix(zPopArg)) zPopArg else matrix(zPopArg, N, zDim, byrow = TRUE)
  list(params = params, zPop = zPop, omega = omega, a = a,
       intercept = intercept, beta = beta, selected = selected,
       covNames = prep$covNames, elboTrace = elboTrace, parHist = parHist,
       mu = last$mu, zPopMat = zPopMat, prep = prep,
       nMix = nMix, mixProb = mixProb, mixnum = last$mixnum)
}

#' Fit entry: prepare data, set up the FOCEi inner problem once, train.
#' @noRd
.vaeFitModel <- function(env) {
  .ui <- env$ui
  .control <- if (exists("vaeControl", envir = env)) env$vaeControl else vaeControl()
  .prep <- .vaeDataPrep(.ui, env$data)
  ## mixture info from the ui: nMix components with probs (p1,...,1-sum)
  .nMix <- tryCatch(as.integer(.ui$saemNMix), error = function(e) 1L)
  if (is.na(.nMix) || .nMix < 1L) .nMix <- 1L
  .mixProb <- 1
  if (.nMix > 1L) {
    .p <- as.numeric(.prep$th[.ui$thetaMixIndex])
    .mixProb <- c(.p, 1 - sum(.p))
  }
  ## parameter-history / iteration-print names: structural typical values on the
  ## mu-referenced etas, the omega diagonal, and the residual error params
  .map <- .foceiEtaThetaMap(.ui)
  .structIdx <- which(!.prep$isFree)
  .parInfo <- list(structIdx = .structIdx,
                   structNames = .map$thetaForEta[.structIdx],
                   omegaNames = paste0("o(", .prep$etaNames, ")"),
                   aNames = names(.prep$a))
  ## back-transform codes for the printed walk (X row: exp/expit/probit thetas)
  .parInfo$xform <- .iterPrintXParFromUi(
    .ui, c(.parInfo$structNames, .parInfo$omegaNames, .parInfo$aNames))
  ## set up the inner likelihood once (compiled model + processed data)
  .innerEnv <- .vaeInnerSetup(.ui, env$data, matrix(0, .prep$N, .prep$zDim), .control)
  on.exit(.vaeInnerFree(), add = TRUE)
  .vaeTrain(.prep, .innerEnv, .control, .nMix, .mixProb, parInfo = .parInfo)
}
