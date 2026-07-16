## Mixture (sub-population) support for NPAG -- a mix() model splits each subject
## into per-component pseudo-subjects; NPAG marginalizes the conditional
## likelihood over the components using the mixture proportions from the mix()
## model: p(y_i | phi) = sum_m mixProb_m * p(y_i | phi, component m).
## Real fit -> weekly slow batch.

nmTest({
  .mixMod <- function() {
    ini({ tka <- log(1.5); tv <- log(32); tcl1 <- log(2); tcl2 <- log(0.5)
      eta.ka ~ 0.3; eta.v ~ 0.1; p1 <- 0.5; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v)
      cl <- mix(exp(tcl1), p1, exp(tcl2)); ke <- cl / v
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - ke * center
      cp <- center / v; cp ~ add(add.sd) })
  }
  .mixMod9 <- function() {
    ini({ tka <- log(1.5); tv <- log(32); tcl1 <- log(2); tcl2 <- log(0.5)
      eta.ka ~ 0.3; eta.v ~ 0.1; p1 <- 0.9; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v)
      cl <- mix(exp(tcl1), p1, exp(tcl2)); ke <- cl / v
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - ke * center
      cp <- center / v; cp ~ add(add.sd) })
  }
  .ctl <- function() npagControl(points = 32L, cycles = 1L, seed = 1L,
                                 gammaOptimize = FALSE)

  test_that("est='npag' fits a mix() model and marginalizes over the components", {
    f <- nlmixr2(.mixMod, nlmixr2data::theo_sd, est = "npag", control = .ctl())
    expect_s3_class(f, "nlmixr2FitData")
    expect_true(is.finite(as.numeric(f$objf)))
  })

  test_that("est='npag' mixture objective depends on the mixture proportion p1", {
    # marginalization is active: changing the mixing proportion must move the
    # nonparametric marginal likelihood (a mixProb-independent objf would mean
    # the components were never combined -- the pre-fix bug).
    f5 <- nlmixr2(.mixMod, nlmixr2data::theo_sd, est = "npag", control = .ctl())
    f9 <- nlmixr2(.mixMod9, nlmixr2data::theo_sd, est = "npag", control = .ctl())
    expect_false(isTRUE(all.equal(as.numeric(f5$objf), as.numeric(f9$objf))))
  })
})
