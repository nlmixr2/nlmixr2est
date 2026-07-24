## Omega off-diagonal (correlated eta) estimation across the newer estimation
## methods.  A correlated 2-eta block `eta.cl + eta.v ~ c(...)` must be
## ESTIMATED (the reported off-diagonal moves off its ini value), not silently
## frozen.  Real fits -> weekly slow batch.

nmTest({
  .omCorMod <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.cl + eta.v ~ c(0.1,
                         0.01, 0.1)
      add.sd <- 0.7 })
    model({ ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) })
  }

  .expectOffDiagEstimated <- function(f) {
    .om <- f$omega
    expect_equal(dim(.om), c(2L, 2L))
    expect_equal(.om[1L, 2L], .om[2L, 1L])
    # the correlation was estimated, not frozen at its 0.01 ini value
    expect_false(isTRUE(all.equal(.om[1L, 2L], 0.01)))
    # and stays a valid covariance
    expect_true(all(is.finite(.om)))
    expect_true(all(eigen(.om, symmetric = TRUE, only.values = TRUE)$values > 0))
  }

  test_that("est='npag' estimates the omega off-diagonal", {
    f <- nlmixr2(.omCorMod, nlmixr2data::theo_sd, est = "npag",
                 control = npagControl(points = 256L, cycles = 15L,
                                       gammaOptimize = FALSE))
    .expectOffDiagEstimated(f)
  })

  test_that("est='npb' estimates the omega off-diagonal", {
    f <- nlmixr2(.omCorMod, nlmixr2data::theo_sd, est = "npb",
                 control = npbControl(points = 50L, burnin = 100L, nsamp = 100L))
    .expectOffDiagEstimated(f)
  })
})
