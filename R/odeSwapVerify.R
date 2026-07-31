# R-visible introspection for the shared solve pool (src/odeSwap.cpp).  These
# exist so tests can assert the MECHANISM -- which model sized the pool, whether
# the private lhs buffer was needed, whether an override leaked into a later fit
# -- rather than only that the numbers matched.

#' Pool decision for the currently registered ODE models.
#'
#' Counters that survive teardown (`pinnedN`, `pinCalledN`, `pooledSolveN`)
#' let a test tell a working pooled solve from a silent rxSolve fallback.
#'
#' @return list with `models` (a data.frame of slot/name/neq/nlhs/loaded/
#'   sizesPool/deny), the chosen `poolSlot`/`poolName`/`poolNeq`/`poolNlhs`,
#'   `maxNlhs`/`maxNlhsSlot`, `scratchNlhs`/`needsScratch`, `overrideNeeded`,
#'   `nLoaded`, and the live `opNeq`/`opNlhs`.
#' @noRd
.odeSwapInfo <- function() odeSwapInfo_()

#' Pool decision for a hypothetical set of models.
#'
#' Pure: takes the per-model state and lhs counts directly, so the tie-break and
#' the "widest lhs is not the largest model" case can be tested without a fit.
#'
#' @param neq integer vector of ODE state counts (0 = slot not loaded)
#' @param nlhs integer vector of lhs output counts, same length as `neq`
#' @return list as in [.odeSwapInfo()], restricted to the plan fields
#' @noRd
.odeSwapPlanFor <- function(neq, nlhs) {
  odeSwapPlanFor_(as.integer(neq), as.integer(nlhs))
}

#' Slot names in registry order, matching the C++ `OdeSwapSlot` enum.
#'
#' The last three are the analytic path's augmented models: the order-2 gradient
#' model, an order-1 model for AGQ quadrature nodes, and a covariance model over
#' its own direction set.  They coexist for the life of a fit -- each has its own
#' compiled entry points, and solving one calls those rather than re-registering
#' a shared slot.
#' @noRd
.odeSwapSlots <- c("inner", "pred", "thetaSens", "hess2",
                   "outer", "outerNode", "outerCov")

#' Drive the shared bad-solve retry loop with stub side effects.
#'
#' Exercises `odeSwapRetryCore` -- the real loop the FOCEi inner, theta-sens,
#' analytic-outer and nlm solves all use -- without needing an ODE that fails.
#' The first `nFail` solves report bad, then they succeed.
#' @noRd
.odeSwapRetryTest <- function(nFail, maxOdeRecalc = 5L, stickyRecalcN = 4L,
                              odeRecalcFactor = 10^0.5, relaxMode = 1L,
                              sticky0 = 0L, restoreTolOnSuccess = TRUE) {
  odeSwapRetryTest_(as.integer(nFail), as.integer(maxOdeRecalc),
                    as.integer(stickyRecalcN), as.double(odeRecalcFactor),
                    as.integer(relaxMode), as.integer(sticky0),
                    isTRUE(restoreTolOnSuccess))
}

#' Relaxation modes, matching the C++ `OdeRelaxMode` enum.
#' @noRd
.odeRelaxGlobal <- 0L
.odeRelaxInd <- 1L

#' Count of proven calc_lhs width mismatches (bound code vs declared nlhs).
#'
#' Non-zero means some model's compiled entry point does not compute the columns the
#' model declares -- see nlmixr2/rxode2#1171.  Cheap enough to read around a single
#' pooled solve.
#' @return integer count since the last registry reset
#' @noRd
.odeSwapLhsMismatchN <- function() {
  .n <- tryCatch(odeSwapInfo_()$lhsWidthMismatchN, error = function(e) NA_real_)
  if (is.null(.n) || length(.n) != 1L || is.na(.n)) return(0L)
  as.integer(.n)
}
