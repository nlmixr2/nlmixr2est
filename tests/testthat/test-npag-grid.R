## Initial grid + condensation for the nonparametric engines.  Core, always-run:
## pure numerics, no model fit.  The Sobol generator is not bit-matched to
## Pmetrics (which uses sobol_burley) -- the converged NPAG solution is robust to
## the low-discrepancy covering -- so the grid is property-tested; the
## deterministic weight-threshold is checked against an in-R reference and the QR
## pass against its independence property.

test_that("npSobolGrid is in-box, deterministic, and covers the box", {
  lo <- c(-2, -1, -3); hi <- c(2, 1, 3)
  g <- npSobolGrid_(2028L, lo, hi)
  expect_equal(dim(g), c(2028L, 3L))
  for (j in 1:3) expect_true(all(g[, j] >= lo[j] & g[, j] <= hi[j]))
  expect_identical(g, npSobolGrid_(2028L, lo, hi))   # deterministic
  # low-discrepancy: coordinate means near the box centre, spread near full range
  expect_true(max(abs(colMeans(g) - (lo + hi) / 2)) < 0.05)
  for (j in 1:3) {
    expect_lt(min(g[, j]), lo[j] + 0.05 * (hi[j] - lo[j]))
    expect_gt(max(g[, j]), hi[j] - 0.05 * (hi[j] - lo[j]))
  }
})

test_that("npCondense weight-threshold matches the in-R reference (Yamada Alg 3)", {
  set.seed(1); psi <- matrix(runif(12 * 10, 0.1, 2), 12, 10)
  lam <- c(0.3, 0.25, 0.2, 0.15, 0.05, 1e-5, 1e-6, 0.02, 1e-7, 0.03)
  ratio <- 1e-3
  ref <- which(lam > max(lam) * ratio)
  r <- npCondense_(lam, psi, ratio, 1e-8)
  expect_equal(r$weightKeep, ref)
})

test_that("npCondense QR drops linearly dependent support points", {
  set.seed(2); base <- matrix(runif(12 * 5, 0.1, 2), 12, 5)
  # columns 6,7 are exact duplicates of 1,2 -> rank 5
  psi <- cbind(base, base[, 1:2])
  lam <- rep(1 / 7, 7)
  r <- npCondense_(lam, psi, 1e-3, 1e-8)
  expect_length(r$qrKeep, 5L)
  expect_equal(qr(psi[, r$qrKeep])$rank, 5L)     # kept set is full rank
  expect_true(all(r$qrKeep %in% seq_len(7)))
})

nmTest({
  test_that(".npEtaBox builds a symmetric auto box from the omega SD", {
    mod <- function() {
      ini({ tka <- log(1.5); tv <- log(32); tke <- log(0.08)
        eta.ka ~ 0.25; eta.v ~ 0.04; eta.ke ~ 0.09
        add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); ke <- exp(tke + eta.ke)
        d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - ke * center
        cp <- center / v; cp ~ add(add.sd) })
    }
    ui <- rxode2::assertRxUi(mod)
    b <- .npEtaBox(ui, npagControl())
    expect_equal(b$names, c("eta.ka", "eta.v", "eta.ke"))
    # auto width = 4 * sd; sd = sqrt(c(0.25,0.04,0.09)) = c(0.5,0.2,0.3)
    expect_equal(b$upper, 4 * c(0.5, 0.2, 0.3), tolerance = 1e-8)
    expect_equal(b$lower, -b$upper, tolerance = 1e-8)
    # grid built over the box is in-box
    g <- npSobolGrid_(200L, b$lower, b$upper)
    for (j in 1:3) expect_true(all(g[, j] >= b$lower[j] & g[, j] <= b$upper[j]))
  })
})
