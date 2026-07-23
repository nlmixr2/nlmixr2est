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

    ## "auto" over the threshold but L0Learn missing: error, not a slow fallback,
    ## and it names the covSelectMaxExact=Inf escape hatch
    expect_error(nlmixr2est:::.vaeCovSelectModes(nCand, ctl("auto"), avail = FALSE),
                 "covSelectMaxExact=Inf")

    ## "auto" below the threshold with L0Learn missing: exact search, no error
    m <- nlmixr2est:::.vaeCovSelectModes(c(5L, 10L), ctl("auto"), avail = FALSE)
    expect_identical(m$mode, c(0L, 0L))
    expect_identical(m$used, "bnb")

    ## covSelectMaxExact=Inf forces the exact search even for a wide problem, so a
    ## missing L0Learn is fine
    m <- nlmixr2est:::.vaeCovSelectModes(nCand, ctl("auto", Inf), avail = FALSE)
    expect_identical(m$mode, c(0L, 0L, 0L))
    expect_identical(m$used, "bnb")
    m <- nlmixr2est:::.vaeCovSelectModes(nCand, ctl("auto", Inf), avail = TRUE)
    expect_identical(m$mode, c(0L, 0L, 0L))  # Inf beats availability

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

  test_that("vaeControl validates covSelectMaxExact", {
    ## Inf and whole numbers pass through
    expect_identical(vaeControl(covSelectMaxExact = Inf)$covSelectMaxExact, Inf)
    expect_identical(vaeControl(covSelectMaxExact = 20L)$covSelectMaxExact, 20L)
    expect_identical(vaeControl(covSelectMaxExact = 20)$covSelectMaxExact, 20L)
    ## a finite non-integer is rejected, not silently truncated to 17
    expect_error(vaeControl(covSelectMaxExact = 17.9))
    ## out-of-range / bad values are rejected
    expect_error(vaeControl(covSelectMaxExact = 0))
    expect_error(vaeControl(covSelectMaxExact = -5))
    expect_error(vaeControl(covSelectMethod = "bogus"))
  })

  test_that("the run-time message stays inside the runInfo one-line budget", {
    ## $runInfo renders one bullet per warning; CLAUDE.md caps these at 75 chars
    ctl <- list(covSelectMethod = "auto", covSelectMaxExact = 25L)
    msgs <- nlmixr2est:::.vaeCovSelectModes(30L, ctl, avail = TRUE)$msg
    expect_gt(length(msgs), 0L)
    expect_true(all(nchar(msgs) <= 75L))
  })

  ## ---- candidate-scoring kernel (vaeScoreSupports_) -------------------------

  ## exact objective, recomputed in-test so the kernel is checked against a
  ## reference that does not call it
  scoreOf <- function(sel, y, X, omega, penalty) {
    Xs <- cbind(1, X[, sel, drop = FALSE])
    coef <- qr.solve(Xs, y)
    sum((y - Xs %*% coef)^2) / omega + penalty * length(sel)
  }
  ## every subset of 0..(nCov-1), 0-based, as vaeScoreSupports_ wants them
  allSupports <- function(nCov) {
    lapply(0:(2^nCov - 1), function(m) which(bitwAnd(m, bitwShiftL(1L, 0:(nCov - 1))) > 0) - 1L)
  }

  test_that("scoring the full enumeration reproduces the exact search", {
    ## the candidate path must differ from the branch-and-bound ONLY in which
    ## subsets it looks at -- same OLS, same score, same tie-break
    .testSeed(21L)
    for (rep in 1:6) {
      N <- 80L; nCov <- sample(4:9, 1L)
      X <- matrix(rnorm(N * nCov), N, nCov)
      k <- sample(0:3, 1L)
      sel <- if (k > 0) sort(sample.int(nCov, k)) else integer(0)
      y <- as.numeric(0.7 + (if (k > 0) X[, sel, drop = FALSE] %*% runif(k, 1, 3) else 0) +
                        rnorm(N, sd = 0.5))
      omega <- 0.4; penalty <- log(N)
      ref <- vaeBestSubset_(matrix(y, ncol = 1), X, omega, FALSE, penalty)
      got <- vaeScoreSupports_(y, X, omega, penalty, allSupports(nCov), polish = FALSE)
      expect_identical(as.integer(got$selected), as.integer(ref$selected[1, ]),
                       info = paste0("rep ", rep))
      expect_equal(got$intercept, as.numeric(ref$intercept[1]), tolerance = 1e-10)
      expect_equal(as.numeric(got$beta), as.numeric(ref$beta[1, ]), tolerance = 1e-10)
    }
  })

  test_that("candidate scoring ignores malformed support indices", {
    .testSeed(22L)
    N <- 60L; nCov <- 5L
    X <- matrix(rnorm(N * nCov), N, nCov)
    y <- as.numeric(X[, 2] * 2 + rnorm(N, sd = 0.3))
    omega <- 0.5; penalty <- log(N)
    clean <- vaeScoreSupports_(y, X, omega, penalty, list(integer(0), 1L), polish = FALSE)
    ## duplicated, out-of-range and negative entries collapse to the same support
    dirty <- vaeScoreSupports_(y, X, omega, penalty,
                               list(integer(0), c(1L, 1L, 99L, -3L)), polish = FALSE)
    expect_identical(dirty$selected, clean$selected)
    expect_equal(dirty$beta, clean$beta, tolerance = 1e-12)
  })

  test_that("the local search never worsens the incumbent and reaches the optimum", {
    .testSeed(23L)
    for (rep in 1:8) {
      N <- 90L; nCov <- sample(5:10, 1L)
      X <- matrix(rnorm(N * nCov), N, nCov)
      k <- sample(1:3, 1L)
      sel <- sort(sample.int(nCov, k))
      y <- as.numeric(0.3 + X[, sel, drop = FALSE] %*% runif(k, 1.5, 3) + rnorm(N, sd = 0.4))
      omega <- 0.5; penalty <- log(N)
      ## deliberately useless candidate set: only the intercept-only model
      bare <- list(integer(0))
      noPolish <- vaeScoreSupports_(y, X, omega, penalty, bare, polish = FALSE)
      polished <- vaeScoreSupports_(y, X, omega, penalty, bare, polish = TRUE)
      sNo <- scoreOf(which(noPolish$selected == 1L), y, X, omega, penalty)
      sYes <- scoreOf(which(polished$selected == 1L), y, X, omega, penalty)
      expect_lte(sYes, sNo)
      ## and from that bare start it still lands on the exact optimum here
      ref <- vaeBestSubset_(matrix(y, ncol = 1), X, omega, FALSE, penalty)
      expect_identical(as.integer(polished$selected), as.integer(ref$selected[1, ]),
                       info = paste0("rep ", rep))
    }
  })

  test_that("L0Learn candidates plus polish match the exact optimum", {
    skip_if_not_installed("L0Learn")
    .testSeed(24L)
    N <- 120L
    for (rep in 1:10) {
      nCov <- sample(8:14, 1L)
      corr <- rep %% 2L == 0L          # half correlated, half independent
      X <- matrix(rnorm(N * nCov), N, nCov)
      if (corr) {                      # rho ~ 0.7 common factor
        f <- rnorm(N)
        X <- sqrt(0.7) * matrix(f, N, nCov) + sqrt(0.3) * X
      }
      k <- sample(1:4, 1L)
      sel <- sort(sample.int(nCov, k))
      y <- as.numeric(0.5 + X[, sel, drop = FALSE] %*% runif(k, 1.5, 3) *
                        sample(c(-1, 1), k, TRUE) + rnorm(N, sd = 0.6))
      omega <- 0.5; penalty <- log(N)
      cand <- nlmixr2est:::.vaeL0Supports(X, y)
      got <- vaeScoreSupports_(y, X, omega, penalty, cand, polish = TRUE)
      ref <- vaeBestSubset_(matrix(y, ncol = 1), X, omega, FALSE, penalty)
      expect_identical(as.integer(got$selected), as.integer(ref$selected[1, ]),
                       info = paste0("rep ", rep, " nCov ", nCov, " corr ", corr))
    }
  })

  test_that(".vaeL0Candidates proposes per dimension in the reduced design", {
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
    ## no pinning: reduced design == full design, so column 4 is index 3
    expect_true(any(vapply(cand[[1]], function(e) identical(e, 3L), logical(1))))

    ## pinCovariates: indices are into the REDUCED design, so they stay in
    ## 0..(length(allowed)-1) -- the C++ caller maps them back to global columns
    allowed <- list(c(3L, 8L), integer(0))
    cand <- nlmixr2est:::.vaeL0Candidates(y, covMat, mode = c(1L, 1L), allowed = allowed)
    expect_true(all(unlist(cand[[1]]) %in% c(0L, 1L)))
    expect_identical(cand[[2]], list(integer(0)))
  })
})
