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

#' Tracked population parameters (structural typical values, omega diagonal,
#' residual error) as a single named vector -- one parameter-history row.
#' @noRd
.vaeParRow <- function(zPop, omega, a, parInfo) {
  .z <- if (length(parInfo$structIdx)) setNames(zPop[parInfo$structIdx], parInfo$structNames) else numeric(0)
  c(.z, setNames(omega, parInfo$omegaNames), setNames(as.numeric(a), parInfo$aNames))
}

#' Encode the per-eta error type of each residual-error parameter for the C++
#' training loop's closed-form error M-step: 0 = additive, 1 = proportional,
#' 2 = other (kept at its current value). Order matches `prep$errType`.
#' @noRd
.vaeErrTypeCode <- function(errType) {
  vapply(errType, function(t) if (identical(t, "add")) 0L else if (identical(t, "prop")) 1L else 2L,
         integer(1), USE.NAMES = FALSE)
}

#' Train the VAE: burn-in (encoder-only, tiny KL) -> main EM (KL anneal + M-step)
#' -> EMA smoothing. The heavy loop runs entirely in C++ (`vaeTrainCpp_`); this
#' function only prepares the inputs (encoder init, prep-derived buffers,
#' iteration-print names/back-transform codes) and re-shapes the return. The
#' inner FOCEi problem must already be set up (`.vaeInnerSetup`) -- the C++ loop
#' drives it through the same likInner0/lpInner engine, per gradient step,
#' without re-running foceiSetup_.
#' @noRd
.vaeTrain <- function(prep, innerEnv, control, nMix = 1L, mixProb = 1,
                      parInfo = NULL) {
  ## RNG is seeded ONCE for the whole estimation in nlmixr2Est.vae (rxWithSeed),
  ## which also covers the model's own random draws and restores the caller's seed
  zDim <- prep$zDim; hDim <- control$hiddenDim; nCov <- ncol(prep$covIn); N <- prep$N
  sigma0 <- if (is.null(control$sigma0)) rep(0.1, zDim) else rep_len(control$sigma0, zDim)
  params <- .vaeEncoderInitParams(zDim, hDim, nCov, prep$zPop, sigma0)

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
  .row0 <- .vaeParRow(prep$zPop, prep$omega, prep$a, parInfo)

  ## prep buffers the C++ loop needs, in the layout vaeTrainCpp_ unpacks: 0-based
  ## theta indices (-1 for a free/mixture eta), error-type codes, per-subject
  ## observed DV, and the plain-matrix covariate design.
  prepC <- c(prep, list(
    zPopThetaIdx0 = ifelse(is.na(prep$zPopThetaIdx), -1L, as.integer(prep$zPopThetaIdx) - 1L),
    errThetaIdx0 = as.integer(prep$errThetaIdx) - 1L,
    errTypeCode = .vaeErrTypeCode(prep$errType),
    yList = lapply(prep$subj, function(s) as.numeric(s$y))))

  .cores <- tryCatch({
    .c <- control$rxControl$cores
    if (is.null(.c) || is.na(.c) || .c < 1L) as.integer(rxode2::getRxThreads()) else as.integer(.c)
  }, error = function(e) 1L)

  ## the print level lives in iterPrintControl$every (absorbed by vaeControl);
  ## surface it as control$print for the C++ loop's final parHist print gate
  control$print <- as.integer(control$iterPrintControl$every)

  .fit <- vaeTrainCpp_(params, prepC, control, as.integer(nMix), as.numeric(mixProb),
                       .cores, .row0, names(.row0), control$iterPrintControl,
                       parInfo$xform, as.integer(parInfo$structIdx) - 1L)

  .selected <- matrix(as.logical(.fit$selected), zDim, ncol(prep$covMat))
  list(params = .fit$params, zPop = as.numeric(.fit$zPop), omega = as.numeric(.fit$omega),
       a = setNames(as.numeric(.fit$a), names(prep$a)),
       intercept = as.numeric(.fit$intercept), beta = .fit$beta, selected = .selected,
       covNames = prep$covNames, elboTrace = as.numeric(.fit$elboTrace), parHist = .fit$parHist,
       mu = .fit$mu, zPopMat = .fit$zPopMat, prep = prep,
       nMix = nMix, mixProb = mixProb, mixnum = as.integer(.fit$mixnum))
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
