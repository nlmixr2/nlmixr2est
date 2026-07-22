## Shared checks for the residual-parameter estimator (plans/vae-els-residual.md).
##
## Written against a specific failure mode: during this work, FOUR wrong
## implementations passed a suite whose assertions were all of the form
## "the objective improved".  Each of them did improve the objective:
##
##   * a joint structural+residual solve that diverged to a corner,
##   * a transform-both-sides objective that transformed `f` as well as `dv`,
##   * a hand-rolled variance that duplicated (and diverged from) the model's own,
##   * a residual that COLLAPSED TO ZERO -- and still beat the closed-form
##     estimator on objective value while doing it.
##
## So an objective comparison is necessary but nowhere near sufficient.  These
## helpers assert the things that actually distinguish a right answer from a
## wrong one: the estimate MOVED, it is not sitting on a bound, and where a
## closed form exists the optimizer reproduces it.

## Every free residual parameter must be strictly interior to its bounds.
## A parameter resting on a bound is the signature of a degenerate optimum --
## a zero variance looks finite to the likelihood (r == 0 becomes r = 1) and so
## looks attractive rather than forbidden.
expectResidInterior <- function(fit, prep, tol = 1e-6) {
  nm <- prep$regressNames
  if (length(nm) == 0L) return(invisible(NULL))
  isErr <- prep$regressErrIdx0 >= 0
  for (i in which(isErr)) {
    v <- unname(fit$theta[[nm[i]]])
    lo <- prep$regressLower[i]; hi <- prep$regressUpper[i]
    testthat::expect_true(is.finite(v),
                          info = paste0(nm[i], " is not finite"))
    if (is.finite(lo)) {
      testthat::expect_gt(v, lo + tol * max(1, abs(lo)))
    }
    if (is.finite(hi)) {
      testthat::expect_lt(v, hi - tol * max(1, abs(hi)))
    }
  }
  invisible(NULL)
}

## The estimate must differ from where it started.  The defect this catches is a
## silent freeze: an unhandled error form was left at its ini() value, which any
## finiteness- or objective-based check passes.
expectMovedFromIni <- function(fit, ui, names, tol = 1e-3) {
  idf <- ui$iniDf
  for (n in names) {
    ini <- idf$est[match(n, idf$name)]
    testthat::expect_gt(abs(unname(fit$theta[[n]]) - ini), tol,
                        label = paste0(n, " moved off its ini() value"))
  }
  invisible(NULL)
}

## Where the optimum is known in closed form the optimizer must find it.  This
## is the assertion that has actually caught objective bugs: an objective that
## is wrong in shape cannot land on sqrt(SSE/n) by accident.
expectMatchesClosedForm <- function(fitOpt, fitMoment, name, tolerance = 5e-3) {
  testthat::expect_equal(unname(fitOpt$theta[[name]]),
                         unname(fitMoment$theta[[name]]),
                         tolerance = tolerance,
                         label = paste0(name, " (optimizer)"),
                         expected.label = paste0(name, " (closed form)"))
  invisible(NULL)
}
