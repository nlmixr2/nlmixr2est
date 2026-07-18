## Exact L0/BIC covariate selection (vaeBestSubset_), the branch-and-bound kernel
## that replaced the exhaustive 2^nCov enumeration in the VAE covariate M-step.
## Validated against a golden fixture whose expected supports are the authoritative
## brute-force optimum of the same objective the reference MIQP (cvxpy + GUROBI)
## solves -- so this needs no Python/GUROBI at run time. See
## tools/vaeGenMiqpFixtures.R. Fast (pure linear algebra, no ODE solve): essential.

nmTest({
  ## Independent brute-force oracle, recomputed in-test so the kernel is checked
  ## against a reference that does not read the fixture.
  bruteBestSubset <- function(y, X, omega, penalty) {
    nCov <- ncol(X)
    best <- list(score = Inf, sel = integer(0))
    for (m in 0:(2^nCov - 1)) {
      sel <- which(bitwAnd(m, bitwShiftL(1L, 0:(nCov - 1))) > 0)
      Xs <- cbind(1, X[, sel, drop = FALSE])
      coef <- qr.solve(Xs, y)
      rss <- sum((y - Xs %*% coef)^2)
      score <- rss / omega + penalty * length(sel)
      if (score < best$score - 1e-9) best <- list(score = score, sel = sel)
    }
    best
  }

  test_that("vaeBestSubset_ matches the golden MIQP-equivalent optimum", {
    g <- readRDS(test_path("baselines", "vae-miqp-golden.rds"))
    for (nm in names(g)) {
      ci <- g[[nm]]; ip <- ci$inputs
      omega <- rep_len(ip$omega, ci$zDim)
      isFree <- rep_len(ip$isFree, ci$zDim)
      got <- vaeBestSubset_(ip$mu, ip$covMat, omega, isFree, ci$penalty)
      ## selected support must match exactly (the core selection claim)
      expect_equal(got$selected, ci$expected$selected,
                   info = paste0("selected mismatch in case ", ci$name))
      if (isTRUE(ci$selectionOnly)) next
      ## coefficients on the chosen support match the OLS/MIQP values
      expect_equal(as.numeric(got$intercept), as.numeric(ci$expected$intercept),
                   tolerance = 1e-6, info = paste0("intercept mismatch: ", ci$name))
      expect_equal(as.numeric(got$beta), as.numeric(ci$expected$beta),
                   tolerance = 1e-6, info = paste0("beta mismatch: ", ci$name))
    }
  })

  test_that("branch-and-bound equals independent brute force on random problems", {
    set.seed(11L)
    for (rep in 1:8) {
      N <- 90L; nCov <- sample(4:11, 1L)
      X <- matrix(rnorm(N * nCov), N, nCov)
      k <- sample(0:3, 1L)
      sel <- if (k > 0) sort(sample.int(nCov, k)) else integer(0)
      y <- 0.7 + (if (k > 0) X[, sel, drop = FALSE] %*% runif(k, 1, 3) * sample(c(-1, 1), k, TRUE) else 0) +
        rnorm(N, sd = 0.5)
      omega <- 0.4; penalty <- log(N)
      ref <- bruteBestSubset(as.numeric(y), X, omega, penalty)
      got <- vaeBestSubset_(matrix(y, ncol = 1), X, omega, FALSE, penalty)
      expect_equal(which(got$selected[1, ] == 1L), ref$sel,
                   info = paste0("rep ", rep, " nCov ", nCov))
    }
  })

  test_that("covSelectAlpha is tunable in vaeControl and defaults to the reference 2", {
    expect_equal(vaeControl()$covSelectAlpha, 2)
    expect_equal(vaeControl(covSelectAlpha = 3)$covSelectAlpha, 3)
    ## round-trips through the list -> vaeControl reconstruction used at dispatch
    ctl <- vaeControl(covSelectAlpha = 1.5)
    expect_equal(do.call(vaeControl, unclass(ctl))$covSelectAlpha, 1.5)
    ## the ramp multiplier must be >= 1 (penalty is only ever inflated early)
    expect_error(vaeControl(covSelectAlpha = 0.5))
  })

  test_that("nCov=32 selects the exact sparse support quickly (regresses 1u<<32 overflow)", {
    g <- readRDS(test_path("baselines", "vae-miqp-golden.rds"))
    ci <- g$E; ip <- ci$inputs
    t0 <- proc.time()[["elapsed"]]
    got <- vaeBestSubset_(ip$mu, ip$covMat, ip$omega, ip$isFree, ci$penalty)
    dt <- proc.time()[["elapsed"]] - t0
    ## previously nothing was ever selected (loop ran only the empty set)
    expect_true(any(got$selected[1, ] == 1L))
    expect_equal(which(got$selected[1, ] == 1L), which(ci$expected$selected[1, ] == 1L))
    expect_lt(dt, 10)
  })
})
