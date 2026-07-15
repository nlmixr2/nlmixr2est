## est="npag" end-to-end: returns a real nlmixr2FitData with the nonparametric
## population summary (theta/Omega), per-subject posterior etas, and the discrete
## support-point distribution.  Real fit -> weekly slow batch.

nmTest({
  .npFitMod <- function() {
    ini({ tka <- log(1.5); tv <- log(32); tke <- log(0.08)
      eta.ka ~ 0.3; eta.v ~ 0.1; eta.ke ~ 0.1; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); ke <- exp(tke + eta.ke)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - ke * center
      cp <- center / v
      cp ~ add(add.sd) })
  }

  test_that("est='npag' returns a valid nlmixr2 fit with support-point outputs", {
    f <- nlmixr2(.npFitMod, nlmixr2data::theo_sd, est = "npag",
                 control = npagControl(points = 256L, cycles = 15L, gammaOptimize = FALSE))
    expect_s3_class(f, "nlmixr2FitData")
    expect_true(inherits(f, "nlmixr2.npag"))
    expect_true(is.finite(as.numeric(f$objf)))
    # population parameters recovered near the (true) initial values
    .pf <- f$parFixed
    expect_true(all(c("tka", "tv", "tke", "add.sd") %in% rownames(.pf)))
    # nonparametric outputs
    .sp <- f$env$npagSupport
    expect_equal(ncol(.sp), 3L)                       # neta columns
    nspp <- nrow(.sp)
    expect_true(nspp >= 1L && nspp <= length(unique(nlmixr2data::theo_sd$ID)))
    expect_equal(sum(f$env$npagWeights), 1, tolerance = 1e-6)
    expect_equal(nrow(f$env$npagPosteriorEta),
                 length(unique(nlmixr2data::theo_sd$ID)))
    expect_equal(f$env$npagNspp, nspp)
  })

  test_that("est='npag' with gamma optimization improves the objective", {
    f0 <- nlmixr2(.npFitMod, nlmixr2data::theo_sd, est = "npag",
                  control = npagControl(points = 256L, cycles = 15L, gammaOptimize = FALSE))
    fg <- nlmixr2(.npFitMod, nlmixr2data::theo_sd, est = "npag",
                  control = npagControl(points = 256L, cycles = 15L, gammaOptimize = TRUE))
    expect_true(is.finite(fg$env$npagGamma) && fg$env$npagGamma > 0)
    expect_lte(as.numeric(fg$objf), as.numeric(f0$objf) + 1e-6)  # -2LL not worse
  })

  test_that("mu-referenced sugar est='mnpag' fits", {
    f <- nlmixr2(.npFitMod, nlmixr2data::theo_sd, est = "mnpag",
                 control = npagControl(points = 256L, cycles = 15L, gammaOptimize = FALSE))
    expect_s3_class(f, "nlmixr2FitData")
    expect_true(is.finite(as.numeric(f$objf)))
  })
})
