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
