# Simulation for the Panhard & Samson (2009) two-level SAEM reproduction.
#
# Panhard X, Samson A (2009).  "Extension of the SAEM algorithm for nonlinear
# mixed effects models with two levels of random effects."  Biostatistics
# 10(1), 121-135.  arXiv:0803.4437.
#
# Section 4.1 design, verbatim where the paper states it:
#   one-compartment first-order absorption,
#     f(t, phi) = D*Ka / (V*Ka - Cl) * (exp(-Cl/V * t) - exp(-Ka * t))
#   phi = (log V, log Ka, log AUC), AUC = D/Cl, D = 4
#   times 0.25, 0.5, 1, 2, 3.5, 5, 7, 9, 12, 24; K = 2 periods; J = 10 per period
#   mu    = (log V = -0.73, log Ka = 0.39, log AUC = 4.61)
#   Omega = diag(0.01, 0.04, 0.04)      (between-subject)
#   Psi   = diag(0.0025, 0.01, 0.01)    (within-subject / inter-occasion)
#   "A combined error model is assumed by setting g(t,phi) = 1 + f(t,phi)"
#   sigma^2 = 0.01
#   n = 24 and n = 40, 1000 replications
#
# UNRESOLVED: the paper reports bias/RMSE for beta_V, beta_Ka and beta_AUC but
# the arXiv HTML never states the beta values used to generate the data.  This
# harness therefore simulates under H0 (beta = 0), which is what the type I
# error study needs and what makes mu/Omega/Psi/sigma comparable; the beta rows
# of the paper's table are out of reach until those values are recovered from
# the published PDF.

panhardTruth <- list(
  mu = c(lV = -0.73, lKa = 0.39, lAUC = 4.61),
  omega = c(lV = 0.01, lKa = 0.04, lAUC = 0.04),
  psi = c(lV = 0.0025, lKa = 0.01, lAUC = 0.01),
  sigma = 0.1,
  dose = 4,
  times = c(0.25, 0.5, 1, 2, 3.5, 5, 7, 9, 12, 24),
  nOcc = 2L
)

#' The paper's structural model
#'
#' @param t times
#' @param v,ka,cl parameters on the natural scale
#' @param dose amount
#' @return concentration
panhardF <- function(t, v, ka, cl, dose = panhardTruth$dose) {
  dose * ka / (v * ka - cl) * (exp(-cl / v * t) - exp(-ka * t))
}

#' Simulate one replicate of the paper's design
#'
#' Each period is a separate dosing episode (`EVID = 4` resets the system), so
#' the two occasions are independent exactly as the crossover design intends.
#'
#' @param n number of subjects
#' @param seed RNG seed for this replicate
#' @param truth list of true values; defaults to the paper's
#' @return data frame with ID, TIME, DV, AMT, EVID and occ
panhardSim <- function(n, seed, truth = panhardTruth) {
  set.seed(seed)
  .k <- truth$nOcc
  .tm <- truth$times
  .j <- length(.tm)
  .out <- vector("list", n * .k)
  .idx <- 0L
  for (.i in seq_len(n)) {
    .b <- stats::rnorm(3L, 0, sqrt(truth$omega))
    for (.occ in seq_len(.k)) {
      .c <- stats::rnorm(3L, 0, sqrt(truth$psi))
      .phi <- truth$mu + .b + .c
      .v <- exp(.phi[["lV"]])
      .ka <- exp(.phi[["lKa"]])
      .cl <- truth$dose / exp(.phi[["lAUC"]])
      .f <- panhardF(.tm, .v, .ka, .cl, truth$dose)
      # g(t, phi) = 1 + f(t, phi); eps ~ N(0, sigma^2)
      .dv <- .f + (1 + .f) * stats::rnorm(.j, 0, truth$sigma)
      .idx <- .idx + 1L
      .out[[.idx]] <- data.frame(
        ID = .i,
        TIME = c(0, .tm) + (.occ - 1L) * 1000,
        DV = c(NA_real_, .dv),
        AMT = c(truth$dose, rep(NA_real_, .j)),
        EVID = c(4L, rep(0L, .j)),
        occ = .occ)
    }
  }
  .d <- do.call(rbind, .out)
  .d[order(.d$ID, .d$TIME), ]
}

#' The model, written with the two-level term the way a user would
# nolint start: object_usage_linter. rxode2 ini()/model() are NSE blocks:
# every assignment here is model specification, consumed by rxode2, not a
# local variable lintr can see used.
panhardModel <- function() {
  ini({
    tlV <- -0.73
    tlKa <- 0.39
    tlAUC <- 4.61
    add.sd <- 0.1
    prop.sd <- 0.1
    eta.lV ~ 0.01
    eta.lKa ~ 0.04
    eta.lAUC ~ 0.04
    iov.lV ~ 0.0025 | occ
    iov.lKa ~ 0.01 | occ
    iov.lAUC ~ 0.01 | occ
  })
  model({
    lV <- tlV + eta.lV + iov.lV
    lKa <- tlKa + eta.lKa + iov.lKa
    lAUC <- tlAUC + eta.lAUC + iov.lAUC
    v <- exp(lV)
    ka <- exp(lKa)
    cl <- 4 / exp(lAUC)
    cp <- linCmt()
    cp ~ combined1(add.sd, prop.sd)
  })
}
# nolint end
