## saem-style mu-expansion for npag: a non-mu-referenced structural fixed-effect
## theta (err==NA, no eta -- e.g. `ke <- exp(tke)`) is not a grid dimension and
## feeds the ODE, so npag cannot estimate it directly.  The mu-expansion injects a
## pseudo-eta (`ke <- exp(tke + eta.tke)`) so the parameter becomes a grid
## dimension estimated by the support-point distribution (same convention as every
## npag theta: the reference theta stays put and the location lives in the
## support).  Real fit -> weekly slow batch.

nmTest({
  # deliberately-wrong ke start (0.30); theo's ke is ~0.08-0.10
  .mod <- function() {
    ini({ tka <- log(1.5); tv <- log(32); tke <- log(0.30)
      eta.ka ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); ke <- exp(tke)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - ke * center
      cp <- center / v; cp ~ add(add.sd) })
  }
  .ctl <- function(muExpand = TRUE)
    npagControl(points = 128L, cycles = 30L, gammaOptimize = FALSE, seed = 1L,
                calcTables = FALSE, muExpand = muExpand)

  test_that("est='npag' mu-expansion makes a non-mu structural theta estimable", {
    f <- nlmixr2(.mod, nlmixr2data::theo_sd, est = "npag", control = .ctl(TRUE))
    expect_true("eta.tke" %in% rownames(f$omega))
    # effective estimate = reference theta + support-point mean (npag convention)
    .sp <- f$env$npagSupport; .wt <- f$env$npagWeights
    .j <- ncol(.sp)                                   # eta.tke is the last eta
    .keEff <- exp(as.numeric(f$theta[["tke"]]) + sum(.wt * .sp[, .j]))
    expect_true(.keEff > 0.05 && .keEff < 0.15)        # recovers theo ke from a 0.30 start
  })

  test_that("est='npag' muExpand=FALSE leaves the non-mu structural theta held", {
    f <- nlmixr2(.mod, nlmixr2data::theo_sd, est = "npag", control = .ctl(FALSE))
    expect_false("eta.tke" %in% rownames(f$omega))
    expect_equal(as.numeric(f$theta[["tke"]]), log(0.30))   # held at ini
  })
})
