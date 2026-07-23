## L0Learn-backed candidate generation for the VAE covariate M-step
## (R/vaeCovSelectL0.R).  L0Learn only PROPOSES supports; the exact objective
## RSS_S/omega + penalty*|S| still picks the winner, so these tests check the
## proposal contract (shape, dedup, determinism, degradation) and the per-eta
## mode resolution.  Pure linear algebra, no ODE solve: essential.

nmTest({

  test_that(".vaeL0Supports proposes well-formed, deduplicated supports", {
    skip_if_not_installed("L0Learn")
    .testSeed(1L)
    N <- 120L; p <- 12L
    X <- matrix(rnorm(N * p), N, p)
    y <- as.numeric(0.5 + X[, c(3, 7, 10)] %*% c(1.5, -2, 1) + rnorm(N, sd = 0.5))
    s <- nlmixr2est:::.vaeL0Supports(X, y)

    expect_type(s, "list")
    expect_gt(length(s), 1L)
    ## every support is an ascending 0-based integer vector inside range
    for (e in s) {
      expect_type(e, "integer")
      expect_false(is.unsorted(e))
      if (length(e)) {
        expect_gte(min(e), 0L)
        expect_lt(max(e), p)
      }
    }
    ## the intercept-only model is always a candidate
    expect_true(any(vapply(s, length, integer(1)) == 0L))
    ## deduplicated
    expect_equal(anyDuplicated(vapply(s, paste, character(1), collapse = ",")), 0L)
    ## the true support is proposed (0-based)
    expect_true(any(vapply(s, function(e) identical(e, c(2L, 6L, 9L)), logical(1))))
  })

  test_that(".vaeL0Supports is deterministic and does not touch the RNG", {
    skip_if_not_installed("L0Learn")
    .testSeed(2L)
    N <- 90L; p <- 8L
    X <- matrix(rnorm(N * p), N, p)
    y <- as.numeric(X[, 2] * 2 + rnorm(N))
    .seedBefore <- .Random.seed
    a <- nlmixr2est:::.vaeL0Supports(X, y)
    expect_identical(.Random.seed, .seedBefore)  # no RNG consumption
    b <- nlmixr2est:::.vaeL0Supports(X, y)
    expect_identical(a, b)
  })

  test_that(".vaeL0Supports degrades to the empty support on unusable input", {
    skip_if_not_installed("L0Learn")
    .empty <- list(integer(0))
    expect_identical(nlmixr2est:::.vaeL0Supports(matrix(0, 10, 0), rnorm(10)), .empty)
    expect_identical(nlmixr2est:::.vaeL0Supports(matrix(1, 2, 3), rnorm(2)), .empty)
    X <- matrix(rnorm(60), 20, 3); X[1, 1] <- NA_real_
    expect_identical(nlmixr2est:::.vaeL0Supports(X, rnorm(20)), .empty)
    X <- matrix(rnorm(60), 20, 3)
    y <- rnorm(20); y[3] <- Inf
    expect_identical(nlmixr2est:::.vaeL0Supports(X, y), .empty)
  })

  test_that(".vaeCovSelectModes resolves the per-eta mode", {
    ctl <- function(method, maxExact = 25L) {
      list(covSelectMethod = method, covSelectMaxExact = maxExact)
    }
    nCand <- c(30L, 10L, 30L)

    ## explicit "bnb": nothing switches, whatever the size
    m <- nlmixr2est:::.vaeCovSelectModes(nCand, ctl("bnb"), avail = TRUE)
    expect_identical(m$mode, c(0L, 0L, 0L))
    expect_identical(m$used, "bnb")
    expect_length(m$msg, 0L)

    ## "auto" + available: only the dimensions at/over the threshold switch
    m <- nlmixr2est:::.vaeCovSelectModes(nCand, ctl("auto"), avail = TRUE)
    expect_identical(m$mode, c(1L, 0L, 1L))
    expect_identical(m$used, "mixed")
    expect_match(m$msg, "not the exact search", all = FALSE)

    ## "auto" below the threshold everywhere: unchanged, silent
    m <- nlmixr2est:::.vaeCovSelectModes(c(5L, 10L), ctl("auto"), avail = TRUE)
    expect_identical(m$mode, c(0L, 0L))
    expect_identical(m$used, "bnb")
    expect_length(m$msg, 0L)

    ## "auto" over the threshold but L0Learn missing: exact search anyway, warned
    m <- nlmixr2est:::.vaeCovSelectModes(nCand, ctl("auto"), avail = FALSE)
    expect_identical(m$mode, c(0L, 0L, 0L))
    expect_identical(m$used, "bnb")
    expect_match(m$msg, "install L0Learn", all = FALSE)

    ## explicit "l0learn": every dimension with candidates switches
    m <- nlmixr2est:::.vaeCovSelectModes(c(3L, 0L), ctl("l0learn"), avail = TRUE)
    expect_identical(m$mode, c(1L, 0L))
    expect_identical(m$used, "mixed")

    ## explicit "l0learn" without the package is an error naming it
    expect_error(nlmixr2est:::.vaeCovSelectModes(nCand, ctl("l0learn"), avail = FALSE),
                 "L0Learn")

    ## the threshold is honored
    m <- nlmixr2est:::.vaeCovSelectModes(c(24L, 25L), ctl("auto"), avail = TRUE)
    expect_identical(m$mode, c(0L, 1L))
  })

  test_that("every run-time message stays inside the runInfo one-line budget", {
    ## $runInfo renders one bullet per warning; CLAUDE.md caps these at 75 chars
    ctl <- list(covSelectMethod = "auto", covSelectMaxExact = 25L)
    msgs <- c(nlmixr2est:::.vaeCovSelectModes(30L, ctl, avail = TRUE)$msg,
              nlmixr2est:::.vaeCovSelectModes(30L, ctl, avail = FALSE)$msg)
    expect_gt(length(msgs), 0L)
    expect_true(all(nchar(msgs) <= 75L))
  })

  test_that(".vaeL0Candidates maps reduced supports back to global columns", {
    skip_if_not_installed("L0Learn")
    .testSeed(3L)
    N <- 100L; p <- 10L
    covMat <- matrix(rnorm(N * p), N, p)
    y <- cbind(as.numeric(covMat[, 4] * 2 + rnorm(N, sd = 0.3)),
               as.numeric(covMat[, 9] * -3 + rnorm(N, sd = 0.3)))

    ## dim 1 uses L0Learn, dim 2 stays on the exact search -> NULL
    cand <- nlmixr2est:::.vaeL0Candidates(y, covMat, mode = c(1L, 0L))
    expect_length(cand, 2L)
    expect_null(cand[[2]])
    expect_true(any(vapply(cand[[1]], function(e) identical(e, 3L), logical(1))))

    ## restricted (pinCovariates) candidate columns: only allowed globals appear
    allowed <- list(c(3L, 8L), integer(0))
    cand <- nlmixr2est:::.vaeL0Candidates(y, covMat, mode = c(1L, 1L), allowed = allowed)
    expect_true(all(unlist(cand[[1]]) %in% c(3L, 8L)))
    expect_identical(cand[[2]], list(integer(0)))
  })
})
