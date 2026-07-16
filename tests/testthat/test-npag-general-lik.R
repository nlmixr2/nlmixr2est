## Generalized (non-normal) likelihood support for NPAG.  The npag Psi sums the
## inner per-observation llikObs, which for a non-normal endpoint is exactly the
## user's log-likelihood -- so the nonparametric objective is already correct.
## The residual/likelihood parameters (iniDf$err non-NA, e.g. a t-distribution's
## df) are estimated with the same frozen-ODE bounded step as the residual params;
## the ODE freeze is valid because those parameters feed only the post-solve f/r,
## never the state derivatives.  npb (Gibbs) still rejects non-normal endpoints.
## Real fit -> weekly slow batch.

nmTest({
  .tMod <- function() {
    ini({ tka <- log(1.5); tv <- log(32); tke <- log(0.08)
      eta.ka ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7; nu <- 8 })
    model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); ke <- exp(tke)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - ke * center
      cp <- center / v
      cp ~ add(add.sd) + dt(nu) })          # Student-t residual: a general likelihood
  }

  test_that("est='npag' fits a generalized (Student-t) likelihood and estimates its params", {
    f <- nlmixr2(.tMod, nlmixr2data::theo_sd, est = "npag",
                 control = npagControl(points = 48L, cycles = 8L, seed = 1L))
    expect_s3_class(f, "nlmixr2FitData")
    expect_true(is.finite(as.numeric(f$objf)))
    # the t-distribution parameters (err-tagged: add / t) are estimated, not held
    expect_true(as.numeric(f$theta[["add.sd"]]) > 0)
    expect_true(as.numeric(f$theta[["nu"]]) > 0)
    expect_false(isTRUE(all.equal(as.numeric(f$theta[["nu"]]), 8)))   # moved off the start
    # gamma is meaningless for a non-normal endpoint (r == 1) -> forced off
    expect_false(isTRUE(f$control$npGammaOptimize))
    # only err-tagged params are frozen-optimized, so the ODE freeze stays valid
    expect_true(isTRUE(f$control$npResidFreeze))
  })

  test_that("est='npag' with muExpand=FALSE holds a non-mu structural param and notes it", {
    # tke has no eta and is not mu-referenced; with the mu-expansion disabled npag
    # holds it at its ini value and surfaces that (not silently) in the fit runInfo.
    f <- nlmixr2(.tMod, nlmixr2data::theo_sd, est = "npag",
                 control = npagControl(points = 24L, cycles = 2L, seed = 1L,
                                       calcTables = FALSE, muExpand = FALSE))
    expect_false("eta.tke" %in% rownames(f$omega))
    expect_true(any(grepl("structural fixed-effect", f$env$runInfo)))
  })

  test_that("est='npag' with muExpand=TRUE grid-estimates the non-mu structural param", {
    f <- nlmixr2(.tMod, nlmixr2data::theo_sd, est = "npag",
                 control = npagControl(points = 24L, cycles = 2L, seed = 1L,
                                       calcTables = FALSE, muExpand = TRUE))
    # opt-in muExpand injects a pseudo-eta so tke becomes a grid dimension
    expect_true("eta.tke" %in% rownames(f$omega))
    expect_false(any(grepl("structural fixed-effect", f$env$runInfo)))
  })

  test_that("est='npb' still rejects a generalized likelihood", {
    expect_error(
      nlmixr2(.tMod, nlmixr2data::theo_sd, est = "npb",
              control = npbControl(points = 10L, burnin = 5L, nsamp = 5L)),
      "does not support generalized")
  })
})
