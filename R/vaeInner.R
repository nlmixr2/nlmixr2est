# vaeInner.R -- drive the FOCEi inner likelihood directly for the VAE. The inner
# problem is set up ONCE (foceiSetup_ via vaeInnerSetup_), then the per-subject
# (and per-mixture-component) objective/gradient are evaluated through
# likInner0/lpInner in a parallel OpenMP loop (vaeInnerLik) -- reusing the inner
# logic (mixtures, multiple endpoints, error structures, log-likelihood,
# censoring) without the nlmixr2 R interface.

#' A foceiControl carrying the VAE's chosen inner likelihood + solving options.
#' focei -> interaction=1; foce/focep -> interaction=0 (focep = FOCE+, R at the
#' live conditional eta); laplace -> the Laplace method.
#' @noRd
.vaeInnerFoceiControl <- function(control) {
  .lik <- control$likelihood
  .interaction <- if (.lik %in% c("foce", "focep")) 0L else 1L
  .foce <- if (identical(.lik, "focep")) "foce+" else "nonmem"
  foceiControl(rxControl = control$rxControl, maxOuterIterations = 0L,
               maxInnerIterations = 0L, covMethod = "", interaction = .interaction,
               foce = .foce,
               sumProd = control$sumProd, optExpression = control$optExpression,
               literalFix = control$literalFix, literalFixRes = control$literalFixRes,
               addProp = control$addProp, calcTables = FALSE, compress = FALSE,
               eventSens = control$eventSens, indTolRelax = control$indTolRelax,
               maxOdeRecalc = control$maxOdeRecalc, odeRecalcFactor = control$odeRecalcFactor,
               stickyRecalcN = control$stickyRecalcN, print = 0L)
}

#' Set up the FOCEi inner problem for `ui` at its current ini() estimates.
#' Builds the focei opt env (which enriches the control with the model-derived
#' neta/ntheta/foceiMuGroup*/... values), augments the fit-flow-derived fields
#' (needOptimHess from the endpoint distribution, est, nF, printTop, AGQ off),
#' preprocesses the data, and calls the C++ vaeInnerSetup_ (foceiSetup_ +
#' updateTheta). Returns the setup env (keep it alive until vaeInnerFree()).
#' @noRd
.vaeInnerSetup <- function(ui, data, etaMat, control, est = "focei") {
  .ui <- rxode2::rxUiDecompress(ui)
  .fc <- .vaeInnerFoceiControl(control)
  .fc$est <- est
  .ui$control <- .fc
  .env <- .ui$foceiOptEnv
  .env$ui <- .ui
  .env$est <- est
  .env$table <- NULL
  .foceiPreProcessData(data, .env, .ui, .fc$rxControl)
  ## fit-flow-derived control fields
  .env$control$est <- est
  .env$control$printTop <- FALSE
  if (is.null(.env$control$nF)) .env$control$nF <- 0L
  .env$control$needOptimHess <- isTRUE(any(.ui$predDfFocei$distribution != "norm"))
  ## AGQ off
  .env$aqn <- 0L; .env$qx <- double(0); .env$qw <- double(0); .env$qfirst <- FALSE
  .env$nAGQ <- 0L; .env$aqLow <- -Inf; .env$aqHi <- Inf; .env$nEstOmega <- 0L
  .env$etaMat <- etaMat
  vaeInnerSetup_(.env)
  .env
}

#' Free the inner-problem state set up by .vaeInnerSetup.
#' @noRd
.vaeInnerFree <- function() invisible(vaeInnerFree_())

#' Evaluate the inner objective (and optionally the eta-gradient) at `etaMat`
#' (rows = ids: nSub, or nSub*nMix for mixtures) through the parallel C++ driver.
#' @noRd
.vaeInnerEval <- function(etaMat, control, grad = FALSE, preds = FALSE) {
  .cores <- tryCatch({
    .c <- control$rxControl$cores
    if (is.null(.c) || is.na(.c) || .c < 1L) as.integer(rxode2::getRxThreads()) else as.integer(.c)
  }, error = function(e) 1L)
  vaeInnerLik(as.matrix(etaMat), .cores, isTRUE(grad), isTRUE(preds))
}

#' Re-set up the inner problem at new population parameters (theta + omega)
#' without recompiling: reuses the env's compiled inner model and processed data,
#' rebuilds rxInv from the new omega, and re-runs foceiSetup_ + updateTheta.
#' @param env the setup env from .vaeInnerSetup
#' @param theta full theta vector (ntheta order): structural intercepts, error,
#'   covariate betas, mixture probs
#' @param omega diagonal random-effect variances (eta order)
#' @param etaMat starting etas [nsub, neta]
#' @noRd
.vaeInnerUpdate <- function(env, theta, omega, etaMat, diagXform = "sqrt") {
  env$thetaIni <- setNames(as.numeric(theta), paste0("THETA[", seq_along(theta), "]"))
  .om <- diag(omega, length(omega))
  .nm <- env$etaNames
  if (!is.null(.nm) && length(.nm) == nrow(.om)) dimnames(.om) <- list(.nm, .nm)
  env$rxInv <- rxode2::rxSymInvCholCreate(mat = .om, diag.xform = diagXform)
  env$etaMat <- etaMat
  vaeInnerSetup_(env)
  invisible(env)
}

#' One ELBO evaluation using the FOCEi inner likelihood (reused, not
#' reimplemented). likInner0(eta) = p(x|z) + p(z); foceiInnerLp = its eta-grad.
#' loss = p_x_z + alphaKL*(p_z - q_z); the encoder upstream is
#' gZ = foceiInnerLp + (alphaKL-1)*Omega^-1 eta, gLogSigma = -alphaKL. Mixtures
#' (nMix>1) evaluate nSub*nMix ids and combine with -2 logsumexp over mixProb.
#' @noRd
.vaeElboStepInner <- function(params, prep, innerEnv, zPop, omega, a, alphaKL, eps,
                              control, nMix = 1L, mixProb = 1, withGrad = TRUE) {
  N <- prep$N; zDim <- prep$zDim
  ## zPop may be a shared vector or an N x zDim matrix of subject-specific KL
  ## centers (covariate model). The inner problem's THETA baseline is the
  ## structural intercept; its built-in eta prior is centered there, so pxz is
  ## recovered by removing that baseline-centered prior, and the KL uses zPopMat.
  if (is.matrix(zPop)) { zPopMat <- zPop; baseline <- colMeans(zPop) }
  else { zPopMat <- matrix(zPop, N, zDim, byrow = TRUE); baseline <- zPop }
  z0 <- matrix(0, N, zDim)
  fw <- vaeEncoderFwdBwd(prep$dataIn, prep$lengths, prep$covIn, eps,
                         params$Wih, params$Whh, params$bih, params$bhh,
                         params$fcW, params$fcB, zDim, z0, z0)
  mu <- fw$mu; logSigma <- fw$logSigma; L <- fw$L; z <- fw$z
  eta <- sweep(z, 2, baseline)                          # z - baseline [N, zDim]

  ## re-setup the inner problem at the current thetas/omega and evaluate
  th <- .vaeBuildTh(prep, baseline, a)
  etaEval <- if (nMix > 1L) do.call(rbind, rep(list(eta), nMix)) else eta
  .vaeInnerUpdate(innerEnv, th, omega, eta)
  r <- .vaeInnerEval(etaEval, control, grad = withGrad, preds = TRUE)

  ln2pi <- log(2 * pi)
  ## inner-problem baseline-centered prior (-loglik), removed to recover pxz
  pzI <- vapply(seq_len(N), function(i) 0.5 * sum(eta[i, ]^2 / omega + log(omega) + ln2pi), numeric(1))
  ## joint = likInner0 = p_x_z + p_z; combine mixtures per physical subject
  if (nMix > 1L) {
    ## per-subject component log-weights: log mixProb[m] - 0.5*obj_im (obj is the
    ## inner -2ll-family joint); a failed solve (NaN/Inf) makes that component
    ## unselectable. Then -2 logsumexp for the value and argmax for the gradient.
    llMat <- vapply(seq_len(N), function(i)
      vapply(seq_len(nMix), function(m) {
        v <- log(mixProb[m]) - 0.5 * r$obj[(m - 1) * N + i]
        if (is.finite(v)) v else -Inf
      }, numeric(1)), numeric(nMix))                    # [nMix, N]
    if (nMix == 1L) llMat <- matrix(llMat, 1L, N)
    joint2 <- vapply(seq_len(N), function(i) {
      ll <- llMat[, i]; mmax <- max(ll)
      if (!is.finite(mmax)) return(NA_real_)
      -2 * (mmax + log(sum(exp(ll - mmax))))
    }, numeric(1))
    jointTot <- sum(joint2, na.rm = TRUE)
    ## gradient and predictions from the best (argmax) component per subject
    best <- vapply(seq_len(N), function(i) {
      w <- which.max(llMat[, i]); if (length(w) == 0L) 1L else w
    }, integer(1))
    lp <- t(vapply(seq_len(N), function(i) r$lp[(best[i] - 1) * N + i, ], numeric(zDim)))
    preds <- lapply(seq_len(N), function(i) r$f[[(best[i] - 1) * N + i]])
  } else {
    jointTot <- sum(r$obj); lp <- r$lp; preds <- r$f
  }
  pxz <- jointTot - sum(pzI)                            # remove baseline prior
  ## VAE KL centered on the (possibly subject-specific) z_pop
  pz <- sum(vapply(seq_len(N), function(i) 0.5 * sum((z[i, ] - zPopMat[i, ])^2 / omega + log(omega) + ln2pi), numeric(1)))
  qz <- sum(vapply(seq_len(N), function(i) 0.5 * sum(eps[i, ]^2 + ln2pi + 2 * logSigma[i, ]), numeric(1)))
  DKL <- pz - qz
  loss <- pxz + alphaKL * DKL

  grads <- NULL
  if (withGrad) {
    ## d(pxz)/dz = lp - eta/omega; KL adds alphaKL*(z - zPopMat)/omega
    gZ <- lp
    for (i in seq_len(N)) gZ[i, ] <- gZ[i, ] - eta[i, ] / omega + alphaKL * (z[i, ] - zPopMat[i, ]) / omega
    gZ[!is.finite(gZ)] <- 0                              # a failed solve skips that subject
    gLS <- matrix(-alphaKL, N, zDim)
    bw <- vaeEncoderFwdBwd(prep$dataIn, prep$lengths, prep$covIn, eps,
                           params$Wih, params$Whh, params$bih, params$bhh,
                           params$fcW, params$fcB, zDim, gZ, gLS)
    grads <- list(Wih = bw$gWih, Whh = bw$gWhh, bih = bw$gbih, bhh = bw$gbhh,
                  fcW = bw$gFcW, fcB = bw$gFcB)
  }
  list(loss = loss, pxz = pxz, DKL = DKL, grads = grads, mu = mu, L = L, z = z,
       preds = preds, mixnum = if (nMix > 1L) best else rep(1L, N))
}
