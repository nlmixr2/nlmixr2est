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
.fsaemImh <- function(map, etaCur, nchain, nsweep = 1L, cores = 1L) {
  .neta <- ncol(map$eta)
  .nsub <- nrow(map$eta)
  # lower-triangular L with Gamma_i = L L' (NA-filled where the proposal failed)
  .cholGamma <- t(vapply(seq_len(.nsub), function(i) {
    .g <- map$gamma[[i]]
    .L <- tryCatch(t(chol(.g)), error = function(e) matrix(NA_real_, .neta, .neta))
    as.numeric(.L)
  }, numeric(.neta*.neta)))
  .nAcc <- integer(.nsub)
  .eta <- etaCur
  for (.s in seq_len(nsweep)) {
    .r <- fsaemImhKernel_(.eta, map$eta, .cholGamma, as.integer(nchain), as.integer(cores))
    .eta <- .r$eta
    .nAcc <- .nAcc + .r$nAcc
  }
  list(eta = .eta, nAcc = .nAcc)
}
