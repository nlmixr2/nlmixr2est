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
#' @param cLoop Run the whole E-M loop in C++ (`TRUE`), avoiding the per-iteration R
#'   round-trip -- so a phase of estimation can be extended (more iterations) without R
#'   overhead.  The eta draw uses rxode2's per-thread threefry engine with a deterministic,
#'   niter-independent per-(iteration, subject) seed, so it is thread-safe, reproducible for
#'   any core count, and a longer run reproduces the exact per-iteration prefix of a shorter
#'   run at the same seed.  Covers the additive/proportional/combined/power/TBS residuals,
#'   single-random-effect covariate regression, structural fixed effects, general
#'   log-likelihood (`ll()`) endpoints (including the box-constrained `likLbfgs` refinement
#'   of bounded likelihood parameters), additive/proportional BLQ censoring (M2/M3/M4),
#'   mode-centered importance sampling (`impInflate`), multiple endpoints, and mixtures.
#'   `FALSE` (default) uses the R-driven loop, which additionally covers multi-endpoint
#'   models with covariates and models with a fix()ed typical value / residual / omega --
#'   cases the C++ loop does not yet handle (it silently falls back to the R loop for those).
#' @param likLbfgs For a general log-likelihood (`ll()`) endpoint, refine the
#'   fixed-effect likelihood parameters each iteration by a box-constrained L-BFGS-B
#'   optimization of the importance-weighted observation log-likelihood (mirrors the
#'   saem/saemix ind.fix10 step), respecting the parameter bounds from the model, rather
#'   than the default single damped-Newton re-solve step.  `TRUE` (default) for `ll()`
#'   models; ignored for standard residual-error models.
#' @param lbfgsLmm,lbfgsFactr,lbfgsPgtol,lbfgsMaxIter L-BFGS-B tuning for the `likLbfgs`
#'   likelihood-parameter refinement (number of corrections, convergence `factr`/`pgtol`,
#'   and max iterations).
#' @param print Iteration-print frequency: display the parameter walk (population estimates
#'   + omega, with the back-transformed row) every `print` iterations (saem/focei/vae style).
#'   `0` (default) captures the parameter history silently.  The walk is *always* saved to the
#'   fit object's parameter history (`fit$parHist` / `fit$parHistStacked`) regardless.  May
#'   also be an `iterPrintControl()` object.
#' @param printNcol,useColor Iteration-print formatting (columns per row, ANSI color); passed
#'   through to `iterPrintControl()`.
#' @param ... Ignored (reserved for future options).
#' @return A list of class `rpemControl`.
#' @export
rpemControl <- function(nGauss = 1000L, nMH = 50000L, mhBurn = 5000L,
                        niter = 50L, collect = 15L, seed = 42L,
                        atol = 1e-8, rtol = 1e-8, cores = 1L,
                        impInflate = 0, cLoop = FALSE,
                        likLbfgs = TRUE, lbfgsLmm = 5L, lbfgsFactr = 1e7,
                        lbfgsPgtol = 0, lbfgsMaxIter = 20L,
                        print = 0L, printNcol = NULL, useColor = NULL, ...) {
  .xtra <- list(...)
  .iterPrintControl <- .absorbIterPrintControl(print = print, printNcol = printNcol,
                                               useColor = useColor,
                                               iterPrintControl = .xtra$iterPrintControl)
  .ret <- list(nGauss = as.integer(nGauss), nMH = as.integer(nMH),
               mhBurn = as.integer(mhBurn), niter = as.integer(niter),
               collect = as.integer(collect), seed = as.integer(seed),
               atol = atol, rtol = rtol, cores = as.integer(cores),
               impInflate = as.numeric(impInflate), cLoop = isTRUE(cLoop),
               likLbfgs = isTRUE(likLbfgs), lbfgsLmm = as.integer(lbfgsLmm),
               lbfgsFactr = as.numeric(lbfgsFactr), lbfgsPgtol = as.numeric(lbfgsPgtol),
               lbfgsMaxIter = as.integer(lbfgsMaxIter),
               print = .iterPrintControl$every,
               iterPrintControl = .iterPrintControl)
  class(.ret) <- "rpemControl"
  .ret
}
