nmTest({
  ## exact objective, mirroring vaeBnbLeaf: RSS/omega + penalty*|S|
  .score <- function(y, X, sel, omega, penalty) {
    .d <- cbind(1, X[, sel, drop = FALSE])
    .r <- stats::lsfit(.d, y, intercept = FALSE)$residuals
    sum(.r^2) / omega + penalty * length(sel)
  }
  ## brute force over the FEASIBLE supports only (at most one column per group)
  .oracleFeasible <- function(y, X, group, omega, penalty) {
    nCov <- ncol(X)
    best <- NULL; bestScore <- Inf
    for (m in 0:(2^nCov - 1)) {
      sel <- which(bitwAnd(m, 2^(seq_len(nCov) - 1L)) > 0L)
      if (anyDuplicated(group[sel])) next
      s <- .score(y, X, sel, omega, penalty)
      if (s < bestScore - 1e-12) { bestScore <- s; best <- sel }
    }
    list(sel = best, score = bestScore)
  }

  test_that("a NULL group reproduces the unconstrained search exactly", {
    set.seed(11)
    N <- 60L; nCov <- 8L
    X <- matrix(rnorm(N * nCov), N, nCov)
    y <- as.numeric(1 + X[, c(2, 5)] %*% c(1.5, -2) + rnorm(N))
    a <- vaeBestSubset_(matrix(y, ncol = 1), X, 0.4, FALSE, log(N))
    b <- vaeBestSubset_(matrix(y, ncol = 1), X, 0.4, FALSE, log(N), "lifo", NULL)
    expect_equal(a, b)
    ## an all-singleton group is the same thing stated explicitly
    c3 <- vaeBestSubset_(matrix(y, ncol = 1), X, 0.4, FALSE, log(N), "lifo",
                         seq_len(nCov))
    expect_equal(a, c3)
  })

  test_that("the constrained search equals brute force over feasible supports", {
    set.seed(12)
    N <- 70L
    for (rep in 1:6) {
      nCov <- 8L
      ## two shapes each for 3 covariates, plus 2 ungrouped indicator columns
      group <- c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 5L)
      X <- matrix(rnorm(N * nCov), N, nCov)
      ## make paired columns genuinely correlated, as two shapes of one covariate are
      X[, 2] <- X[, 1] * 0.9 + rnorm(N, sd = 0.4)
      X[, 4] <- X[, 3] * 0.9 + rnorm(N, sd = 0.4)
      y <- as.numeric(0.5 + X[, c(2, 5, 7)] %*% c(1.5, -2, 1) + rnorm(N))
      omega <- 0.4; penalty <- log(N)
      got <- vaeBestSubset_(matrix(y, ncol = 1), X, omega, FALSE, penalty,
                            "lifo", group)
      ref <- .oracleFeasible(y, X, group, omega, penalty)
      expect_equal(which(got$selected[1, ] == 1L), ref$sel, info = rep)
      ## and the constraint actually binds: never two columns of one group
      expect_false(anyDuplicated(group[which(got$selected[1, ] == 1L)]) > 0L)
    }
  })

  test_that("all frontier strategies find the same constrained optimum", {
    set.seed(13)
    N <- 60L; nCov <- 6L
    group <- c(1L, 1L, 2L, 2L, 3L, 3L)
    X <- matrix(rnorm(N * nCov), N, nCov)
    X[, 2] <- X[, 1] * 0.95 + rnorm(N, sd = 0.3)
    y <- as.numeric(1 + X[, c(1, 4)] %*% c(2, -1.5) + rnorm(N))
    ref <- vaeBestSubset_(matrix(y, ncol = 1), X, 0.4, FALSE, log(N), "lifo", group)
    for (st in c("fifo", "lc")) {
      got <- vaeBestSubset_(matrix(y, ncol = 1), X, 0.4, FALSE, log(N), st, group)
      expect_equal(got, ref, info = st)
    }
  })

  test_that("a duplicated column loses to its twin only via the group tie-break", {
    ## exactly collinear twins: identical RSS, so the constraint must keep one
    set.seed(14)
    N <- 50L
    x <- rnorm(N)
    X <- cbind(x, x, rnorm(N))
    group <- c(1L, 1L, 2L)
    y <- as.numeric(1 + 2 * x + rnorm(N, sd = 0.2))
    got <- vaeBestSubset_(matrix(y, ncol = 1), X, 0.4, FALSE, log(N), "lifo", group)
    expect_equal(sum(got$selected[1, 1:2]), 1L)
  })

  test_that("candidate scoring repairs infeasible proposals", {
    set.seed(15)
    N <- 60L
    x1 <- rnorm(N); x2 <- rnorm(N)
    ## column 2 is the weaker shape of covariate 1
    X <- cbind(x1, x1 * 0.5 + rnorm(N, sd = 1.5), x2)
    group <- c(1L, 1L, 2L)
    y <- as.numeric(1 + 2 * x1 - 1.5 * x2 + rnorm(N, sd = 0.3))
    ## propose an infeasible support naming BOTH shapes of covariate 1
    got <- vaeScoreSupports_(y, X, 0.4, log(N), list(c(0L, 1L, 2L)),
                             polish = FALSE, group = group)
    sel <- which(got$selected == 1L)
    expect_false(anyDuplicated(group[sel]) > 0L)
    ## the repair keeps the univariately stronger of the two shapes
    expect_true(1L %in% sel)
    ## and with polish it reaches the exact constrained optimum
    pol <- vaeScoreSupports_(y, X, 0.4, log(N), list(c(0L, 1L, 2L)),
                             polish = TRUE, group = group)
    ref <- vaeBestSubset_(matrix(y, ncol = 1), X, 0.4, FALSE, log(N), "lifo", group)
    expect_equal(which(pol$selected == 1L), which(ref$selected[1, ] == 1L))
  })

  test_that("a malformed group vector is rejected", {
    set.seed(16)
    N <- 40L
    X <- matrix(rnorm(N * 4L), N, 4L)
    y <- rnorm(N)
    expect_error(
      vaeBestSubset_(matrix(y, ncol = 1), X, 0.4, FALSE, log(N), "lifo", 1:3),
      "one entry per covariate column")
  })
})
