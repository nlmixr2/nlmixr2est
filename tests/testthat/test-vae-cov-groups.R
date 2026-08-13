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
    ## and specifically the FIRST twin: on an exact score tie the incumbent is
    ## the lexicographically smaller support, so this is deterministic
    expect_equal(as.integer(got$selected[1, 1:2]), c(1L, 0L))
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
    expect_error(
      vaeBestSubset_(matrix(y, ncol = 1), X, 0.4, FALSE, log(N), "lifo", NULL, 1:3),
      "one entry per covariate column")
  })

  ## ---- atomic blocks (the hockey stick's two arms) ---------------------------

  ## WT's linear column beside its two hockey arms.  The arms sum to the linear
  ## column, so this is exactly the collinear case blocks exist to resolve.
  .hockeyDesign <- function(N, seed) {
    set.seed(seed)
    wt <- runif(N, 40, 140)
    ctr <- stats::median(wt)
    cbind(lin = wt - ctr,
          armLow = (wt < ctr) * (wt - ctr),
          armHi = (wt >= ctr) * (wt - ctr))
  }

  test_that("singleton blocks reproduce the unconstrained search exactly", {
    ## the zero-drift guarantee: saying "every column is its own block" must not
    ## change the tree, the branch order or the tie-break for anyone
    set.seed(31)
    for (rep in 1:8) {
      N <- 70L; nCov <- sample(4:9, 1L)
      X <- matrix(rnorm(N * nCov), N, nCov)
      k <- sample(0:3, 1L)
      sel <- if (k > 0) sort(sample.int(nCov, k)) else integer(0)
      y <- as.numeric(0.7 +
                        (if (k > 0) X[, sel, drop = FALSE] %*% runif(k, 1, 3) else 0) +
                        rnorm(N, sd = 0.5))
      grp <- sample(rep(seq_len(ceiling(nCov / 2)), each = 2), nCov)
      for (g in list(NULL, grp)) {
        a <- vaeBestSubset_(matrix(y, ncol = 1), X, 0.4, FALSE, log(N), "lifo", g)
        b <- vaeBestSubset_(matrix(y, ncol = 1), X, 0.4, FALSE, log(N), "lifo", g,
                            seq_len(nCov))
        expect_equal(a, b, info = paste0("rep ", rep))
      }
    }
  })

  test_that("a block is selected whole, or not at all", {
    ## The motivating case.  span{1, armLow, armHi} == span{1, lin, armLow} at the
    ## same column count and so the same penalty -- an exact tie.  Without blocks
    ## the tie-break lands on a form that does not read as a hockey stick.
    N <- 80L
    X <- .hockeyDesign(N, 5)
    set.seed(5)
    y <- as.numeric(1 + X[, "armLow"] * (-0.03) + X[, "armHi"] * 0.06 +
                      rnorm(N, sd = 0.2))
    grp <- c(1L, 1L, 1L)
    ## unconstrained: lin + one arm, the tie-equivalent parameterization
    free <- vaeBestSubset_(matrix(y, ncol = 1), X, 0.25, FALSE, log(N), "lifo")
    expect_equal(which(free$selected[1, ] == 1L), c(1L, 2L))
    ## group only: at most ONE column, so two slopes are unreachable
    gOnly <- vaeBestSubset_(matrix(y, ncol = 1), X, 0.25, FALSE, log(N), "lifo", grp)
    expect_length(which(gOnly$selected[1, ] == 1L), 1L)
    ## group + block: both arms, never one, never lin beside an arm
    blk <- vaeBestSubset_(matrix(y, ncol = 1), X, 0.25, FALSE, log(N), "lifo",
                          grp, c(1L, 2L, 2L))
    expect_equal(which(blk$selected[1, ] == 1L), c(2L, 3L))
  })

  test_that("the blocked search equals brute force over feasible supports", {
    ## exactness against an oracle that never calls the search: enumerate every
    ## support that is block-complete AND takes at most one block per group
    .oracleBlk <- function(y, X, group, block, omega, penalty) {
      nCov <- ncol(X)
      best <- integer(0); bestScore <- Inf
      for (m in 0:(2^nCov - 1)) {
        sel <- which(bitwAnd(m, 2^(seq_len(nCov) - 1L)) > 0L)
        if (any(vapply(unique(block[sel]),
                       function(b) !all(which(block == b) %in% sel),
                       logical(1)))) next
        bl <- unique(block[sel])
        if (anyDuplicated(vapply(bl, function(b) group[which(block == b)[1L]],
                                 integer(1)))) next
        s <- .score(y, X, sel, omega, penalty)
        if (s < bestScore - 1e-12) { bestScore <- s; best <- sel }
      }
      best
    }
    N <- 80L
    group <- c(1L, 1L, 1L, 2L, 2L, 3L)
    block <- c(1L, 2L, 2L, 3L, 4L, 5L)
    for (rep in 1:8) {
      X <- cbind(.hockeyDesign(N, 200 + rep),
                 log(runif(N, 20, 80) / 50), runif(N, 20, 80) - 50,
                 rbinom(N, 1, 0.4))
      ## rotate the truth so hockey wins, lin wins, a plain covariate wins and
      ## nothing wins -- the last two are the no-false-positive cases
      y <- switch(1L + (rep %% 4L),
                  as.numeric(1 + X[, 2] * (-0.02) + X[, 3] * 0.05 + rnorm(N, sd = 0.3)),
                  as.numeric(1 + X[, 1] * 0.02 + rnorm(N, sd = 0.3)),
                  as.numeric(1 + X[, 6] * 0.8 + rnorm(N, sd = 0.3)),
                  as.numeric(1 + rnorm(N, sd = 0.3)))
      got <- vaeBestSubset_(matrix(y, ncol = 1), X, 0.25, FALSE, log(N), "lifo",
                            group, block)
      expect_equal(which(got$selected[1, ] == 1L),
                   .oracleBlk(y, X, group, block, 0.25, log(N)),
                   info = paste0("rep ", rep))
    }
  })

  test_that("the frontier discipline does not change a blocked selection", {
    ## the prune stays admissible with blocks, so the search is still exact and
    ## the strategy only reorders visitation
    N <- 60L
    X <- cbind(.hockeyDesign(N, 3), matrix(rnorm(N * 2L), N, 2L))
    set.seed(3)
    y <- as.numeric(1 + X[, 2] * (-0.03) + X[, 3] * 0.06 + rnorm(N, sd = 0.2))
    g <- c(1L, 1L, 1L, 2L, 3L); b <- c(1L, 2L, 2L, 3L, 4L)
    r <- lapply(c("lifo", "fifo", "lc"), function(s) {
      vaeBestSubset_(matrix(y, ncol = 1), X, 0.25, FALSE, log(N), s, g, b)
    })
    expect_equal(r[[1]], r[[2]])
    expect_equal(r[[1]], r[[3]])
  })

  test_that("a half-block proposal is completed, not discarded", {
    ## L0Learn knows nothing about blocks, so it can propose one arm.  Dropping
    ## it would make hockey unreachable on the approximate path entirely.
    N <- 80L
    X <- .hockeyDesign(N, 5)
    set.seed(5)
    y <- as.numeric(1 + X[, 2] * (-0.03) + X[, 3] * 0.06 + rnorm(N, sd = 0.2))
    g <- c(1L, 1L, 1L); b <- c(1L, 2L, 2L)
    half <- list(integer(0), 1L)             # 0-based: the low arm alone
    got <- vaeScoreSupports_(y, X, 0.25, log(N), half, polish = FALSE, g, b)
    expect_equal(which(got$selected == 1L), c(2L, 3L))
    ## without blocks the same proposal is never completed: the pair is simply
    ## not reachable, so the lone arm is scored on its own (and here loses to the
    ## intercept-only model)
    plain <- vaeScoreSupports_(y, X, 0.25, log(N), half, polish = FALSE, g)
    expect_lte(length(which(plain$selected == 1L)), 1L)
    ## the polish moves whole blocks, so it reaches the exact optimum from empty
    pol <- vaeScoreSupports_(y, X, 0.25, log(N), list(integer(0)), polish = TRUE, g, b)
    ref <- vaeBestSubset_(matrix(y, ncol = 1), X, 0.25, FALSE, log(N), "lifo", g, b)
    expect_equal(which(pol$selected == 1L), which(ref$selected[1, ] == 1L))
    expect_equal(which(pol$selected == 1L), c(2L, 3L))
  })

  test_that("a group collision is repaired by dropping whole blocks", {
    ## a proposal naming BOTH the linear column and a hockey arm doubles up on
    ## the covariate's group.  The repair has to resolve it by keeping one whole
    ## RELATIONSHIP -- splitting a block would leave half a hockey stick, which
    ## the leaf then rejects, silently losing the candidate.
    N <- 80L
    X <- .hockeyDesign(N, 5)
    set.seed(5)
    y <- as.numeric(1 + X[, 2] * (-0.03) + X[, 3] * 0.06 + rnorm(N, sd = 0.2))
    g <- c(1L, 1L, 1L); b <- c(1L, 2L, 2L)
    got <- vaeScoreSupports_(y, X, 0.25, log(N), list(c(0L, 1L)), polish = FALSE,
                             g, b)
    sel <- which(got$selected == 1L)
    ## whatever survives is a whole block, never a mixture of the two
    expect_true(identical(sel, 1L) || identical(sel, c(2L, 3L)) ||
                  identical(sel, integer(0)))
  })

  test_that("an NA block id makes a column its own block", {
    ## NA is the "no constraint" marker on the R side, matching group
    N <- 60L
    X <- .hockeyDesign(N, 3)
    set.seed(3)
    y <- as.numeric(1 + X[, 1] * 0.03 + rnorm(N, sd = 0.2))
    free <- vaeBestSubset_(matrix(y, ncol = 1), X, 0.25, FALSE, log(N), "lifo")
    naBlk <- vaeBestSubset_(matrix(y, ncol = 1), X, 0.25, FALSE, log(N), "lifo",
                            NULL, c(NA_integer_, NA_integer_, NA_integer_))
    expect_equal(free, naBlk)
  })
})
