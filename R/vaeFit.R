# vaeFit.R -- VAE training: encoder init, one ELBO gradient step (encoder C++
# forward/backward + rxode2 decoder), Adam, closed-form M-step, and the
# burn-in -> KL-anneal -> smoothing schedule. Orchestration is thin R; the heavy
# work is the C++ encoder BPTT (vaeEncoderFwdBwd) and the rxode2 decoder solve.

#' Initialize encoder parameters
#' @noRd
.vaeEncoderInitParams <- function(zDim, hDim, nCov, zPop, sigma0, seed = 1L) {
  set.seed(seed)
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

#' Assemble the full theta vector from current z_pop (structural) + a (residual)
#' @noRd
.vaeBuildTh <- function(prep, zPop, a) {
  th <- prep$th
  th[prep$zPopThetaIdx] <- zPop
  if (!is.na(prep$aThetaIdx)) th[prep$aThetaIdx] <- a
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
  th <- .vaeBuildTh(prep, zPop, a)
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
    eta <- z[i, ] - zPop
    E <- .vaeDecoderSolveSubject(am, th, eta, s$ev, s$times)
    if (is.null(E)) return(NULL)
    px <- .vaeDecoderPxz(E, s$y)
    pxz <- pxz + px$pxz
    gZdec[i, ] <- px$gEta
    preds[[i]] <- E$f
  }

  ## KL: p_z - q_z (z_pop shared across subjects in Milestone A)
  ln2pi <- log(2 * pi)
  pz <- 0; qz <- 0
  for (i in seq_len(N)) {
    pz <- pz + 0.5 * sum((z[i, ] - zPop)^2 / omega + log(omega) + ln2pi)
    qz <- qz + 0.5 * sum(eps[i, ]^2 + ln2pi + 2 * logSigma[i, ])
  }
  DKL <- pz - qz
  loss <- pxz + alphaKL * DKL

  grads <- NULL
  if (withGrad) {
    gZ <- gZdec
    for (i in seq_len(N)) gZ[i, ] <- gZ[i, ] + alphaKL * (z[i, ] - zPop) / omega
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
  ssr <- 0
  for (i in seq_len(N)) ssr <- ssr + sum((prep$subj[[i]]$y - preds[[i]])^2)
  aCur <- sqrt(ssr / prep$Nobs)
  list(zPop = zPop + gamma * (zPopCur - zPop),
       omega = omega + gamma * (omegaCur - omega),
       a = a + gamma * (aCur - a))
}

#' Train the VAE: burn-in (encoder-only, tiny KL) -> main EM (KL anneal + M-step)
#' -> EMA smoothing. Returns the trained encoder + population estimates + traces.
#' @noRd
.vaeTrain <- function(prep, am, control, verbose = FALSE) {
  zDim <- prep$zDim; hDim <- control$hiddenDim; nCov <- ncol(prep$covIn); N <- prep$N
  sigma0 <- if (is.null(control$sigma0)) rep(0.1, zDim) else rep_len(control$sigma0, zDim)
  params <- .vaeEncoderInitParams(zDim, hDim, nCov, prep$zPop, sigma0, seed = control$seed)
  adam <- .vaeAdamInit(params)
  zPop <- prep$zPop; omega <- prep$omega; a <- prep$a
  Lg <- control$nGradStep
  set.seed(control$seed)
  .t <- 0L; last <- NULL

  ## burn-in: encoder-only training with a tiny fixed KL weight
  for (it in seq_len(control$itersBurnIn)) {
    for (l in seq_len(Lg)) {
      eps <- matrix(stats::rnorm(N * zDim), N, zDim)
      st <- .vaeElboStep(params, prep, am, zPop, omega, a, 0.001, eps)
      if (is.null(st)) stop("est=\"vae\" decoder solve failed during burn-in", call. = FALSE)
      .t <- .t + 1L
      .ad <- .vaeAdamStep(params, st$grads, adam, control$burnInLearningRate, .t)
      params <- .ad$params; adam <- .ad$state
      .ms <- .vaeMStep(st$mu, st$L, st$preds, prep, zPop, omega, a, 1)
      zPop <- .ms$zPop; omega <- .ms$omega; a <- .ms$a
      last <- st
    }
  }

  ## main EM
  nMain <- control$iters
  elboTrace <- numeric(nMain)
  zPopTrace <- matrix(0, nMain, zDim); omegaTrace <- matrix(0, nMain, zDim); aTrace <- numeric(nMain)
  for (it in seq_len(nMain)) {
    gamma <- if (it <= control$gammaIter) 1 else 1 / (1 + it - control$gammaIter)
    .ms <- .vaeMStep(last$mu, last$L, last$preds, prep, zPop, omega, a, gamma)
    zPop <- .ms$zPop; omega <- .ms$omega; a <- .ms$a
    alphaKL <- if (it <= control$klWarmup) 0.01 + 0.99 * (it - 1) / max(1, control$klWarmup - 1) else 1
    elbos <- numeric(Lg)
    for (l in seq_len(Lg)) {
      eps <- matrix(stats::rnorm(N * zDim), N, zDim)
      st <- .vaeElboStep(params, prep, am, zPop, omega, a, alphaKL, eps)
      if (is.null(st)) stop("est=\"vae\" decoder solve failed during training", call. = FALSE)
      .t <- .t + 1L
      .ad <- .vaeAdamStep(params, st$grads, adam, control$learningRate, .t)
      params <- .ad$params; adam <- .ad$state
      elbos[l] <- st$pxz + st$DKL
      last <- st
    }
    elboTrace[it] <- mean(elbos); zPopTrace[it, ] <- zPop; omegaTrace[it, ] <- omega; aTrace[it] <- a
    if (verbose && it %% 25 == 0)
      message(sprintf("iter %d/%d  -ELBO=%.2f  zPop=(%s)  a=%.3f", it, nMain,
                      elboTrace[it], paste(round(zPop, 3), collapse = ","), a))
  }

  list(params = params, zPop = zPop, omega = omega, a = a,
       elboTrace = elboTrace, zPopTrace = zPopTrace, omegaTrace = omegaTrace, aTrace = aTrace,
       prep = prep, am = am)
}

#' Fit entry: prepare data, build decoder model, train.
#' @noRd
.vaeFitModel <- function(env) {
  .ui <- env$ui
  .control <- if (exists("vaeControl", envir = env)) env$vaeControl else vaeControl()
  .prep <- .vaeDataPrep(.ui, env$data)
  .am <- .vaeDecoderModel(.ui)
  if (is.null(.am)) stop("est=\"vae\" could not build the decoder sensitivity model", call. = FALSE)
  .vaeTrain(.prep, .am, .control)
}
