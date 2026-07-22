## The covariate M-step solves an L0 best-subset problem that the reference
## implementation (Rohleff et al. 2025) poses as a MIQP and hands to GUROBI:
##
##   min over S:  RSS_S/omega + penalty*|S|
##
## We solve it with a self-contained branch-and-bound instead (no commercial
## solver).  The claim in the documentation is that this returns the SAME
## optimum the MIQP would, not an approximation -- so it has to be checked
## against exhaustive enumeration, not merely against itself.

nmTest({
  ## exhaustive minimizer over all 2^p supports
  .brute <- function(y, X, omega, penalty) {
    p <- ncol(X) - 1L
    best <- list(sel = integer(0), score = Inf)
    for (m in 0:(2^p - 1)) {
      sel <- which(bitwAnd(m, 2^(0:(p - 1))) > 0)
      fit <- .lm.fit(X[, c(1L, 1L + sel), drop = FALSE], y)
      sc <- sum(fit$residuals^2) / omega + penalty * length(sel)
      if (sc < best$score - 1e-12) best <- list(sel = sel, score = sc)
    }
    best
  }

  test_that("the branch-and-bound attains the exact L0 optimum", {
    skip_on_cran()
    set.seed(42)
    for (trial in 1:40) {
      N <- sample(30:200, 1)
      p <- sample(1:7, 1)
      X <- cbind(1, matrix(stats::rnorm(N * p), N, p))
      ## a sparse truth, so the optimum is a genuine subset rather than all-in
      b <- c(stats::rnorm(1), ifelse(stats::runif(p) < 0.4, stats::rnorm(p, 0, 2), 0))
      y <- as.vector(X %*% b + stats::rnorm(N, 0, stats::runif(1, 0.2, 2)))
      omega <- stats::runif(1, 0.001, 2)
      penalty <- stats::runif(1, 0.1, 30)
      bf <- .brute(y, X, omega, penalty)
      ## the solver is exact, so every frontier discipline must agree with it
      for (strat in c("lifo", "fifo", "lc")) {
        r <- vaeBestSubsetL0_(y, X, omega, penalty, strat)
        expect_equal(r$score, bf$score, tolerance = 1e-8,
                     info = paste0("trial ", trial, " strategy ", strat))
      }
    }
  })

  test_that("an all-zero-signal problem selects nothing", {
    set.seed(7)
    X <- cbind(1, matrix(stats::rnorm(200 * 4), 200, 4))
    y <- stats::rnorm(200)
    ## a large penalty against pure noise: the intercept-only model must win
    r <- vaeBestSubsetL0_(y, X, 1, 1e6, "lifo")
    expect_length(r$sel, 0L)
  })
})
