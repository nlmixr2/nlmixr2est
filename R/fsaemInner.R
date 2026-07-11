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
  if (!all(.pred$distribution == "norm")) return(FALSE)         # continuous normal (Hessian path TODO)
  .err <- ui$iniDf$err
  .err <- .err[!is.na(.err)]
  if (!all(.err == "add")) return(FALSE)                        # additive residual only for now
  .cov <- ui$muRefCovariateDataFrame                            # no mu-ref covariates (mprior would vary by subject)
  if (!is.null(.cov) && nrow(.cov) > 0L) return(FALSE)
  if (length(ui$mixProbs) > 0L) return(FALSE)                   # no mixtures yet
  TRUE
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
  .env <- .fsaemInnerSetup(ui, data, matrix(0, .N, .neta), .fc)
  cfg$fsaemInnerEnv <- .env
  cfg$fsaemStep <- function(theta, omega, etaCur, nchain, nsweep = 5L) {
    .fsaemInnerUpdate(.env, theta, omega, matrix(0, .N, .neta))
    .map <- .fsaemInnerMap(.fc, .neta)
    .imh <- .fsaemImh(.map, etaCur, as.integer(nchain), as.integer(nsweep))
    .imh$eta
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
