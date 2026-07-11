# fsaemInner.R -- f-SAEM (Karimi, Lavielle & Moulines 2020) proposal builder.
# Reuses the FOCEi inner problem (foceiSetup_ via vaeInnerSetup_) to optimize the
# per-subject conditional MAP of the random effects and read back the FOCEi inner
# information matrix H = Gamma_i^-1 (the IMH proposal precision).  Unlike the VAE
# path (which only evaluates the inner at supplied etas, maxInnerIterations=0),
# this sets a real inner optimizer so innerOpt1() finds the MAP.

#' foceiControl for the f-SAEM inner MAP + proposal covariance.
#'
#' `fastCov`/`fastLik` map to the inner interaction/foce settings:
#'   - "jacobian": interaction=0 -- H = J' Sigma^-1 J + Omega^-1 (paper Eq 17/19)
#'   - "hessian":  interaction=1 -- Gauss-Newton/Laplace Hessian (paper Eq 13);
#'                 non-normal endpoints force the exact finite-diff inner Hessian.
#' @noRd
.fsaemInnerFoceiControl <- function(control) {
  .cov <- control$fastCov
  .lik <- control$fastLik
  .interaction <- if (identical(.cov, "jacobian")) 0L else 1L
  .foce <- if (identical(.lik, "focep")) "foce+" else "nonmem"
  .maxInner <- control$fastInnerIt
  if (is.null(.maxInner) || is.na(.maxInner) || .maxInner < 1L) .maxInner <- 100L
  foceiControl(rxControl = control$rxControl, maxOuterIterations = 0L,
               maxInnerIterations = as.integer(.maxInner), covMethod = "",
               interaction = .interaction, foce = .foce,
               sumProd = control$sumProd, optExpression = control$optExpression,
               literalFix = control$literalFix,
               addProp = control$addProp, calcTables = FALSE, compress = FALSE,
               eventSens = control$eventSens, indTolRelax = control$indTolRelax,
               maxOdeRecalc = control$maxOdeRecalc, odeRecalcFactor = control$odeRecalcFactor,
               print = 0L)
}

#' Set up the FOCEi inner problem for the f-SAEM proposal at `ui`'s current
#' ini() estimates.  Mirrors `.vaeInnerSetup` but keeps a live inner optimizer.
#' Returns the setup env (keep alive until `.fsaemInnerFree()`).
#' @noRd
.fsaemInnerSetup <- function(ui, data, etaMat, control) {
  .ui <- rxode2::rxUiDecompress(ui)
  .fc <- .fsaemInnerFoceiControl(control)
  .fc$est <- "focei"
  .ui$control <- .fc
  .env <- .ui$foceiOptEnv
  .env$ui <- .ui
  .env$est <- "focei"
  .env$table <- NULL
  .foceiPreProcessData(data, .env, .ui, .fc$rxControl)
  .env$control$est <- "focei"
  .env$control$printTop <- FALSE
  if (is.null(.env$control$nF)) .env$control$nF <- 0L
  # non-normal endpoints require the exact (finite-diff) inner Hessian; a Hessian
  # proposal also wants it even for normal data
  .env$control$needOptimHess <-
    isTRUE(any(.ui$predDfFocei$distribution != "norm")) ||
    identical(control$fastCov, "hessian")
  .env$aqn <- 0L; .env$qx <- double(0); .env$qw <- double(0); .env$qfirst <- FALSE
  .env$nAGQ <- 0L; .env$aqLow <- -Inf; .env$aqHi <- Inf; .env$nEstOmega <- 0L
  .env$etaMat <- etaMat
  vaeInnerSetup_(.env)
  .env
}

#' Re-set up the inner problem at new population parameters (theta + omega)
#' without recompiling.  Same mechanism as `.vaeInnerUpdate`.
#' @noRd
.fsaemInnerUpdate <- function(env, theta, omega, etaMat, diagXform = "sqrt") {
  env$thetaIni <- setNames(as.numeric(theta), paste0("THETA[", seq_along(theta), "]"))
  .om <- diag(omega, length(omega))
  .nm <- env$etaNames
  if (!is.null(.nm) && length(.nm) == nrow(.om)) dimnames(.om) <- list(.nm, .nm)
  env$rxInv <- rxode2::rxSymInvCholCreate(mat = .om, diag.xform = diagXform)
  env$etaMat <- etaMat
  vaeInnerSetup_(env)
  invisible(env)
}

#' Free the f-SAEM inner-problem state.
#' @noRd
.fsaemInnerFree <- function() invisible(vaeInnerFree_())

#' Optimize the per-subject MAP and return the proposal (mean + covariance).
#'
#' @return list with `eta` (nsub x neta MAP), `gamma` (nsub-list of neta x neta
#'   proposal covariance = solve(H)), `hess` (nsub-list of H = Gamma^-1), and
#'   `ok` (per-subject convergence flag).
#' @noRd
.fsaemInnerMap <- function(control, neta) {
  .cores <- tryCatch({
    .c <- control$rxControl$cores
    if (is.null(.c) || is.na(.c) || .c < 1L) as.integer(rxode2::getRxThreads()) else as.integer(.c)
  }, error = function(e) 1L)
  .r <- fsaemInnerMap_(.cores)
  .nsub <- nrow(.r$eta)
  .hess <- lapply(seq_len(.nsub), function(i) matrix(.r$hess[i, ], neta, neta))
  .gamma <- lapply(seq_len(.nsub), function(i) {
    if (.r$ok[i] == 0L) return(matrix(NA_real_, neta, neta))
    .H <- .hess[[i]]
    .g <- try(solve(.H), silent = TRUE)
    if (inherits(.g, "try-error")) matrix(NA_real_, neta, neta) else .g
  })
  list(eta = .r$eta, gamma = .gamma, hess = .hess, ok = .r$ok)
}

#' Is `ui` within the current f-SAEM fast-kernel support envelope?
#'
#' The C++ side reconstructs the inner's theta as [population phi (constant across
#' subjects), single additive residual]; models outside that envelope would get a
#' wrong parameterization, so they degrade to standard SAEM instead.
#' @noRd
.fsaemSupported <- function(ui) {
  .pred <- ui$predDf
  if (is.null(.pred) || length(.pred$cond) != 1L) return(FALSE) # single endpoint
  if (length(ui$mixProbs) > 0L) return(FALSE)                   # no mixtures yet
  # General log-likelihood endpoint (ll() ~ expr, distribution=="LL"): the inner
  # supplies the observation likelihood, so the fast kernel handles it even
  # though plain saem cannot.  It must run throughout (distribution=4 path) --
  # forced in .fsaemGeneralLik / the control.
  if (.pred$distribution == "LL") return(TRUE)
  if (!all(.pred$distribution == "norm")) return(FALSE)         # else continuous normal
  .err <- ui$iniDf$err
  .err <- .err[!is.na(.err)]
  if (!all(.err %in% c("add", "prop"))) return(FALSE)           # additive/proportional/combined residual
  # mu-ref covariates are supported: non-time-varying ones are absorbed into the
  # per-subject mprior data, time-varying ones are kept as inner regressor betas
  # refreshed from the live Plambda each iteration.
  TRUE
}

#' TRUE when the fit uses a general log-likelihood endpoint (distribution=4 path)
#' @noRd
.fsaemGeneralLik <- function(ui) {
  .pred <- ui$predDf
  !is.null(.pred) && length(.pred$cond) == 1L && .pred$distribution == "LL"
}

#' Build the f-SAEM fast-simulation step for a SAEM fit and attach it to `cfg`.
#'
#' Sets up the FOCEi inner once and returns a closure the C++ SAEM loop calls
#' each fast iteration with the current estimate.  The closure re-parameterizes
#' the inner (theta = [structural fixed effects, residual] in THETA order; omega
#' = current diagonal), optimizes the per-subject MAP + proposal covariance, and
#' runs the independent Metropolis-Hastings kernel over the chains.
#' @return `cfg` with `$fsaemStep` (closure) and `$fsaemInnerEnv` (keep-alive).
#' @noRd
#' Phi1 (random-effect) parameter bounds in inner-eta order, from the mu-ref
#' theta iniDf bounds (the same bounds L-BFGS-B uses for phi0).  Used to keep the
#' IMH proposals for a general log-likelihood inside plausible values.
#' @return list(lower, upper, nbd) with L-BFGS-B nbd codes (0 none/1 lo/2 both/3 hi)
#' @noRd
.fsaemPhi1Bounds <- function(ui) {
  .iniDf <- ui$iniDf
  .etaRows <- .iniDf[which(.iniDf$neta1 == .iniDf$neta2), ]
  .etaRows <- .etaRows[order(.etaRows$neta1), ]
  .mr <- ui$muRefDataFrame
  .th <- .mr$theta[match(.etaRows$name, .mr$eta)]
  .lo <- .iniDf$lower[match(.th, .iniDf$name)]
  .hi <- .iniDf$upper[match(.th, .iniDf$name)]
  list(lower = ifelse(is.finite(.lo), .lo, 0),
       upper = ifelse(is.finite(.hi), .hi, 0),
       nbd = as.integer(ifelse(is.finite(.lo) & is.finite(.hi), 2L,
                        ifelse(is.finite(.lo), 1L,
                        ifelse(is.finite(.hi), 3L, 0L)))))
}

.fsaemInstallStep <- function(ui, data, rxControl, cfg) {
  if (!.fsaemSupported(ui)) {
    .minfo(paste0("f-SAEM (fast=TRUE) fast kernel not yet supported for this model ",
                  "(needs a single additive-error continuous endpoint, no covariates/mixtures); ",
                  "running standard SAEM"))
    return(cfg)
  }
  .iniDf <- ui$iniDf
  .neta <- sum(.iniDf$neta1 == .iniDf$neta2, na.rm = TRUE)
  .N <- length(unique(data[[if ("ID" %in% names(data)) "ID" else "id"]]))
  .fc <- list(rxControl = rxControl,
              fastCov = rxode2::rxGetControl(ui, "fastCov", "auto"),
              fastLik = rxode2::rxGetControl(ui, "fastLik", "focei"),
              fastInnerIt = 100L,
              sumProd = rxode2::rxGetControl(ui, "sumProd", FALSE),
              optExpression = rxode2::rxGetControl(ui, "optExpression", TRUE),
              literalFix = rxode2::rxGetControl(ui, "literalFix", FALSE),
              addProp = rxode2::rxGetControl(ui, "addProp", "combined2"),
              eventSens = rxode2::rxGetControl(ui, "eventSens", "jump"),
              indTolRelax = rxode2::rxGetControl(ui, "indTolRelax", TRUE),
              maxOdeRecalc = rxode2::rxGetControl(ui, "maxOdeRecalc", 5L),
              odeRecalcFactor = rxode2::rxGetControl(ui, "odeRecalcFactor", 10^0.5))
  # "auto": Jacobian for continuous single-endpoint normal data, else Hessian
  if (identical(.fc$fastCov, "auto")) {
    .fc$fastCov <- if (all(ui$predDf$distribution == "norm")) "jacobian" else "hessian"
  }
  # bound-respecting reproducible IMH proposals (threefry seeded by the SAEM seed).
  # Boundary clamping is ONLY for general log-likelihood models: for normal
  # (add/prop/combined) data the residual error optimizer does the clamping, so
  # the IMH proposal stays unconstrained there.
  .bounds <- if (.fsaemGeneralLik(ui)) .fsaemPhi1Bounds(ui) else NULL
  .seed <- as.integer(rxode2::rxGetControl(ui, "seed", 99))
  .nRetry <- as.integer(rxode2::rxGetControl(ui, "nRetry", 10L))
  .hasCov <- !is.null(ui$muRefCovariateDataFrame) && nrow(ui$muRefCovariateDataFrame) > 0L
  if (.hasCov) {
    # Covariate path: the time-invariant covariate effect is absorbed into the
    # per-subject prior mean, so the inner is built on the mprior-as-data model
    # (each mu-ref intercept is a per-subject nlmixrMprior* data column, etas
    # kept).  The C++ loop passes the full per-subject mprior_phi1 each iteration.
    .iniPhi <- vapply(ui$muRefDataFrame$theta, function(th)
      ui$iniDf$est[match(th, ui$iniDf$name)], numeric(1))
    .mpri0 <- matrix(.iniPhi, .N, .neta, byrow = TRUE)
    .setup <- .fsaemInnerSetupCov(ui, data, .mpri0, .fc)
    cfg$fsaemInnerEnv <- .setup$env
    cfg$fsaemStep <- function(mpriorMat, ares, bres, omega, plambda, etaCur, nchain, kiter, nsweep = 5L) {
      .fsaemInnerUpdateCov(.setup, mpriorMat, ares, bres, plambda, omega)
      .map <- .fsaemInnerMap(.fc, .neta)
      .imh <- .fsaemImh(.map, etaCur, as.integer(nchain), as.integer(nsweep),
                        mprior = mpriorMat, bounds = .bounds, seed = .seed,
                        nRetry = .nRetry, kiter = as.integer(kiter))
      .imh$eta
    }
    return(cfg)
  }
  # No-covariate path: the inner uses the ui model directly with a constant
  # population phi (mprior is the same for every subject).
  .env <- .fsaemInnerSetup(ui, data, matrix(0, .N, .neta), .fc)
  cfg$fsaemInnerEnv <- .env
  # Inner THETA is in UI ntheta order: structural (mu-referenced) positions take
  # the current population phi (mprior row 1); residual positions take the
  # current additive (ares) / proportional (bres) estimate for their endpoint.
  .thetaDf <- ui$iniDf[!is.na(ui$iniDf$ntheta), c("ntheta", "err", "condition")]
  .thetaDf <- .thetaDf[order(.thetaDf$ntheta), ]
  .nTheta <- nrow(.thetaDf)
  .structPos <- which(is.na(.thetaDf$err))
  .residPos <- which(!is.na(.thetaDf$err))
  .residIsAdd <- .thetaDf$err[.residPos] == "add"
  .residEp <- match(.thetaDf$condition[.residPos], ui$predDf$cond) # 1-based endpoint
  # bound vectors for the C++ orchestration (empty = unbounded)
  .lower <- if (is.null(.bounds)) numeric(0) else as.numeric(.bounds$lower)
  .upper <- if (is.null(.bounds)) numeric(0) else as.numeric(.bounds$upper)
  .nbd   <- if (is.null(.bounds)) integer(0) else as.integer(.bounds$nbd)
  cfg$fsaemStep <- function(mpriorMat, ares, bres, omega, plambda, etaCur, nchain, kiter, nsweep = 5L) {
    .theta <- numeric(.nTheta)
    .theta[.structPos] <- mpriorMat[1, ]
    if (length(.residPos)) {
      .theta[.residPos] <- ifelse(.residIsAdd, ares[.residEp], bres[.residEp])
    }
    .cores <- as.integer(rxode2::getRxThreads())
    if (is.na(.cores) || .cores < 1L) .cores <- 1L
    # whole per-iteration orchestration in C++ (inner update + MAP + IMH)
    fsaemStepCpp_(.env, as.numeric(.theta), as.numeric(omega), mpriorMat, etaCur,
                  as.integer(nchain), as.integer(nsweep), .cores,
                  .lower, .upper, .nbd, as.double(.seed), as.integer(.nRetry),
                  as.integer(kiter))
  }
  cfg
}

#' Run `nchain` independent Metropolis-Hastings sweeps of the f-SAEM kernel.
#'
#' @param map result of `.fsaemInnerMap()` (proposal mean + covariance)
#' @param etaCur chain-major current state, ((nchain*nsub) x neta), row = c*nsub + id
#' @param nchain number of chains
#' @param nsweep number of IMH sweeps to run (each sweep = one proposal per
#'   subject per chain); the SAEM simulation step uses a small number
#' @return list with updated `eta` and per-subject acceptance count `nAcc`
#'   (summed over sweeps and chains)
#' @noRd
.fsaemImh <- function(map, etaCur, nchain, nsweep = 1L, cores = 1L,
                      mprior = NULL, bounds = NULL, seed = 0L, nRetry = 10L,
                      kiter = 0L) {
  .neta <- ncol(map$eta)
  .nsub <- nrow(map$eta)
  # lower-triangular L with Gamma_i = L L' (NA-filled where the proposal failed)
  .cholGamma <- t(vapply(seq_len(.nsub), function(i) {
    .g <- map$gamma[[i]]
    .L <- tryCatch(t(chol(.g)), error = function(e) matrix(NA_real_, .neta, .neta))
    as.numeric(.L)
  }, numeric(.neta*.neta)))
  .mp <- if (is.null(mprior)) matrix(0, .nsub, .neta) else mprior
  .lower <- if (is.null(bounds)) numeric(0) else as.numeric(bounds$lower)
  .upper <- if (is.null(bounds)) numeric(0) else as.numeric(bounds$upper)
  .nbd <- if (is.null(bounds)) integer(0) else as.integer(bounds$nbd)
  .nAcc <- integer(.nsub)
  .eta <- etaCur
  # the kernel is serial but must size the threefry engine array for the full
  # thread count so a later fit's parallel solves see a correctly-sized array
  .cores <- as.integer(rxode2::getRxThreads())
  if (is.na(.cores) || .cores < 1L) .cores <- 1L
  for (.s in seq_len(nsweep)) {
    # Iteration-indexed threefry stream base = f(kiter, sweep), mirroring
    # est="imp" (seed0 + iter*nExp); the kernel adds c*nsub + id per chain/
    # subject.  Being a pure function of (kiter, sweep, chain, subject) -- not a
    # running counter -- makes every iteration independently reproducible, so the
    # number of steps in a phase can be changed dynamically and still reproduced.
    .streamBase <- as.double(seed) +
      (as.double(kiter) * nsweep + (.s - 1)) * (as.double(nchain) * .nsub)
    .r <- fsaemImhKernel_(.eta, map$eta, .cholGamma, as.integer(nchain), .cores,
                          .mp, .lower, .upper, .nbd, .streamBase, as.integer(nRetry))
    .eta <- .r$eta
    .nAcc <- .nAcc + .r$nAcc
  }
  list(eta = .eta, nAcc = .nAcc)
}
