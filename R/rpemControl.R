#' RPEM estimation control
#'
#' Options for the Randomized Parametric Expectation Maximization (RPEM)
#' estimation method (Chen et al. 2024).  This is the K=1 minimal core; see
#' `design/rpem/` for the full roadmap.
#'
#' @param nGauss Number of Monte Carlo samples per subject in the E-step
#'   (`m_Gauss`).
#' @param nMH Number of Metropolis-Hastings trials collected in the M-step.
#' @param mhBurn Number of M-step MH burn-in trials (discarded).
#' @param niter Maximum number of E-M iterations.
#' @param collect Number of terminal iterations averaged for the final estimate.
#' @param seed RNG seed for the threefry sampler (reproducible for a fixed thread
#'   count).
#' @param atol,rtol ODE solver tolerances.
#' @param cores Number of cores for the threefry draw (solve threading is set by
#'   rxode2).
#' @param impInflate Opt-in mode-centered importance sampling for the E-step.  `0`
#'   (default) keeps the paper's prior sampling (draw eta ~ N(0, Omega)).  A value
#'   `>= 1` draws instead from N(EBE, impInflate*Omega) -- centered at the previous
#'   iteration's posterior mean with that variance-inflation factor -- and
#'   importance-weights, improving posterior-tail coverage for high-variance random
#'   effects in multi-eta models (whose largest Omega prior sampling under-estimates).
#'   Experimental: a partial mitigation, not a full fix (see design/rpem/04).
#' @param ... Ignored (reserved for future options).
#' @return A list of class `rpemControl`.
#' @export
rpemControl <- function(nGauss = 1000L, nMH = 50000L, mhBurn = 5000L,
                        niter = 50L, collect = 15L, seed = 42L,
                        atol = 1e-8, rtol = 1e-8, cores = 1L,
                        impInflate = 0, ...) {
  .ret <- list(nGauss = as.integer(nGauss), nMH = as.integer(nMH),
               mhBurn = as.integer(mhBurn), niter = as.integer(niter),
               collect = as.integer(collect), seed = as.integer(seed),
               atol = atol, rtol = rtol, cores = as.integer(cores),
               impInflate = as.numeric(impInflate))
  class(.ret) <- "rpemControl"
  .ret
}
