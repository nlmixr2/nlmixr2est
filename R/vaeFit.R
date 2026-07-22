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
  ## The whole step (encoder fwd/bwd, per-subject decoder solve loop + p(x|z), KL)
  ## runs in C++ (vaeDecoderElboStep_).  Only the augmented-model rxode2 solve stays
  ## in R -- the C++ calls this per-subject closure, which builds the full theta
  ## (baseline structural TVs + residual error) and solves subject i at its etas.
  ## The decoder baseline (THETA) is immaterial to f (only THETA+ETA=z matters), so
  ## use the intercept.
  baseline <- if (is.matrix(zPop)) colMeans(zPop) else zPop
  th <- .vaeBuildTh(prep, baseline, a)
  .etav <- am$dirs
  .solve <- function(i0, e, t) {
    s <- prep$subj[[i0 + 1L]]
    .foceiAnalyticSolveFA(am, c(th, setNames(e, .etav)), s$ev, s$times, tol = t)
  }
  .yList <- lapply(prep$subj, function(s) as.numeric(s$y))
  vaeDecoderElboStep_(params, prep, zPop, as.numeric(omega), as.numeric(a),
                      as.numeric(alphaKL), as.matrix(eps), .solve, .yList,
                      isTRUE(withGrad), 1e-10, 5L, 10^(0.5), TRUE)
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
#' residual error, and any nonMuTheta="regress" fixed-effect thetas) as a single
#' named vector -- one parameter-history row.  `parInfo` supplies only the column
#' metadata (indices/names); the printed structural typical values and omega
#' diagonal skip fixed latent dims (`structIdx`/`omegaIdx`), and the regressed
#' thetas are appended last (values via `regressVals`) so the column order matches
#' the C++ `parRow` lambda and the `xform` codes.
#' @noRd
.vaeParRow <- function(zPop, omega, a, parInfo, regressVals = NULL) {
  .z <- if (length(parInfo$structIdx)) setNames(zPop[parInfo$structIdx], parInfo$structNames) else numeric(0)
  .om <- if (length(parInfo$omegaIdx)) setNames(omega[parInfo$omegaIdx], parInfo$omegaNames) else numeric(0)
  .reg <- if (length(parInfo$regressNames)) setNames(as.numeric(regressVals), parInfo$regressNames) else numeric(0)
  c(.z, .om, setNames(as.numeric(a), parInfo$aNames), .reg)
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
    .sIdx <- which(!prep$isFree & !prep$zPopFix)
    .oIdx <- which(!prep$zPopFix)
    parInfo <- list(structIdx = .sIdx, structNames = prep$etaNames[.sIdx],
                    omegaIdx = .oIdx, omegaNames = paste0("o(", prep$etaNames[.oIdx], ")"),
                    aNames = names(prep$a))
  }
  ## nonMuTheta="regress": surface the regressed fixed-effect thetas in the
  ## parameter-history walk.  parInfo carries only their names (metadata); their
  ## starting VALUES (the ini() estimates at the regressed theta indices, 0-based
  ## `regressThetaIdx0` into the full theta) are passed to .vaeParRow explicitly.
  parInfo$regressNames <- prep$regressNames
  .regressVals0 <- if (length(prep$regressThetaIdx0)) {
    prep$th[prep$regressThetaIdx0 + 1L]
  } else numeric(0)
  .row0 <- .vaeParRow(prep$zPop, prep$omega, prep$a, parInfo, regressVals = .regressVals0)

  ## prep buffers the C++ loop needs, in the layout vaeTrainCpp_ unpacks: 0-based
  ## theta indices (-1 for a free/mixture eta), error-type codes, per-subject
  ## observed DV, and the plain-matrix covariate design.
  prepC <- c(prep, list(
    zPopThetaIdx0 = ifelse(is.na(prep$zPopThetaIdx), -1L, as.integer(prep$zPopThetaIdx) - 1L),
    errThetaIdx0 = as.integer(prep$errThetaIdx) - 1L,
    errTypeCode = .vaeErrTypeCode(prep$errType),
    yList = lapply(prep$subj, function(s) as.numeric(s$y))))
  ## nonMuTheta="regress": 0-based full-theta indices + ini bounds of the fixed
  ## thetas the C++ M-step regresses with bobyqa (empty when not in regress mode)
  prepC$regressThetaIdx0 <- as.integer(prep$regressThetaIdx0)
  prepC$regressLower <- as.numeric(prep$regressLower)
  prepC$regressUpper <- as.numeric(prep$regressUpper)
  ## latent dims whose structural theta is fixed (held at ini by the M-step)
  prepC$zPopFix <- as.logical(prep$zPopFix)
  ## pinned covariate selection: per-(eta x covariate) allow-mask restricting the
  ## branch-and-bound to model-declared pairs.  Drop the NULL placeholder when
  ## pinning is inactive so the C++ containsElementNamed guard sees no mask and
  ## runs the full search.
  prepC$covAllow <- NULL
  if (!is.null(prep$covAllow)) {
    prepC$covAllow <- matrix(as.integer(prep$covAllow), prep$zDim, ncol(prep$covMat))
  }

  .cores <- tryCatch({
    .c <- control$rxControl$cores
    if (is.null(.c) || is.na(.c) || .c < 1L) as.integer(rxode2::getRxThreads()) else as.integer(.c)
  }, error = function(e) 1L)

  ## surface the parallel-encoder-backward non-reproducibility in $runInfo (this
  ## warning is collected into the fit's run information); only relevant when it
  ## actually parallelizes (cores > 1)
  if (isTRUE(control$parEncoderBackward) && .cores > 1L) {
    warning("encoder: small parallel deviation; parEncoderBackward=FALSE turns off",
            call. = FALSE)
  }

  ## the print level lives in iterPrintControl$every (absorbed by vaeControl);
  ## surface it as control$print for the C++ loop's final parHist print gate
  control$print <- as.integer(control$iterPrintControl$every)

  ## nonMuTheta="grad": stash the per-fit context the analytic outer-gradient
  ## M-step reads; the C++ loop then passes only theta/eta/omega per M-step
  if (identical(control$nonMuTheta, "grad") && length(prep$regressNames)) {
    .vaeGradInit(innerEnv$ui, innerEnv$dataSav, prep$regressNames)
    ## .vaeGradEnv lives for the SESSION; drop the fit-specific state when this
    ## fit ends so a later focei fast fit cannot see it (see .foceiAnalyticSolveAll)
    on.exit(.vaeGradReset(), add = TRUE)
  }

  .fit <- vaeTrainCpp_(params, prepC, control, as.integer(nMix), as.numeric(mixProb),
                       .cores, .row0, names(.row0), control$iterPrintControl,
                       parInfo$xform, as.integer(parInfo$structIdx) - 1L)

  .selected <- matrix(as.logical(.fit$selected), zDim, ncol(prep$covMat))
  list(params = .fit$params, zPop = as.numeric(.fit$zPop), omega = as.numeric(.fit$omega),
       a = setNames(as.numeric(.fit$a), names(prep$a)),
       intercept = as.numeric(.fit$intercept), beta = .fit$beta, selected = .selected,
       covNames = prep$covNames, elboTrace = as.numeric(.fit$elboTrace), parHist = .fit$parHist,
       mu = .fit$mu, zPopMat = .fit$zPopMat, prep = prep,
       regressTheta = setNames(as.numeric(.fit$regressTheta), prep$regressNames),
       nRegGrad = as.integer(.fit$nRegGrad), nRegFallback = as.integer(.fit$nRegFallback),
       nMix = nMix, mixProb = mixProb, mixnum = as.integer(.fit$mixnum))
}

#' Fit entry: prepare data, set up the FOCEi inner problem once, train.
#' @noRd
.vaeFitModel <- function(env) {
  .ui <- env$ui
  .control <- if (exists("vaeControl", envir = env)) env$vaeControl else vaeControl()
  .prep <- .vaeDataPrep(.ui, env$data, .control)
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
  ## printed structural + omega columns: the estimated latent-space parameters,
  ## i.e. the dims whose backing theta is NOT fixed.  A fixed theta (e.g.
  ## nonMuTheta="fix") is held at ini with a fixed omega, so BOTH its typical value
  ## (structIdx also drops free/mixture etas) and its omega (omegaIdx) are excluded
  ## from the iteration print.
  .structIdx <- which(!.prep$isFree & !.prep$zPopFix)
  .omegaIdx <- which(!.prep$zPopFix)
  .parInfo <- list(structIdx = .structIdx,
                   structNames = .map$thetaForEta[.structIdx],
                   omegaIdx = .omegaIdx,
                   omegaNames = paste0("o(", .prep$etaNames[.omegaIdx], ")"),
                   aNames = names(.prep$a))
  ## back-transform codes for the printed walk (X row: exp/expit/probit thetas).
  ## Include the nonMuTheta="regress" thetas so their column gets the right
  ## back-transform (they are appended last, matching .vaeParRow / the C++ parRow).
  .parInfo$xform <- .iterPrintXParFromUi(
    .ui, c(.parInfo$structNames, .parInfo$omegaNames, .parInfo$aNames, .prep$regressNames))
  ## set up the inner likelihood once (compiled model + processed data)
  .innerEnv <- .vaeInnerSetup(.ui, env$data, matrix(0, .prep$N, .prep$zDim), .control)
  on.exit(.vaeInnerFree(), add = TRUE)
  .vaeTrain(.prep, .innerEnv, .control, .nMix, .mixProb, parInfo = .parInfo)
}
