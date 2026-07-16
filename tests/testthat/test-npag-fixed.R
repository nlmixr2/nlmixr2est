## est="npag" fixed-parameter support: fixed thetas (incl. residual) are held at
## their ini value; fixed-Omega etas (e.g. IOV) stay support dimensions but keep
## their variance held at the fixed value instead of being estimated.  Real fits
## -> weekly slow batch.

nmTest({

  test_that("est='npag' holds a fixed population theta at its ini value", {
    .mod <- function() {
      ini({ tka <- log(1.5); tv <- fix(log(31.5)); tke <- log(0.08); add.sd <- 0.7
        eta.ka ~ 0.3; eta.ke ~ 0.1 })
      model({ ka <- exp(tka + eta.ka); v <- exp(tv); ke <- exp(tke + eta.ke)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - ke * center
        cp <- center / v; cp ~ add(add.sd) })
    }
    f <- nlmixr2(.mod, nlmixr2data::theo_sd, est = "npag",
                 control = npagControl(points = 200L, cycles = 6L))
    expect_s3_class(f, "nlmixr2FitData")
    # the fixed theta stays exactly at its ini value (log scale)
    expect_equal(unname(f$parFixedDf["tv", "Estimate"]), log(31.5), tolerance = 1e-8)
  })

  test_that("est='npag' holds a fixed residual (add.sd) at its ini value", {
    .mod <- function() {
      ini({ tka <- log(1.5); tv <- log(31.5); tke <- log(0.08); add.sd <- fix(0.7)
        eta.ka ~ 0.3; eta.ke ~ 0.1 })
      model({ ka <- exp(tka + eta.ka); v <- exp(tv); ke <- exp(tke + eta.ke)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - ke * center
        cp <- center / v; cp ~ add(add.sd) })
    }
    f <- nlmixr2(.mod, nlmixr2data::theo_sd, est = "npag",
                 control = npagControl(points = 200L, cycles = 6L))
    expect_s3_class(f, "nlmixr2FitData")
    expect_equal(unname(f$parFixedDf["add.sd", "Estimate"]), 0.7, tolerance = 1e-8)
  })

  test_that("est='npag' fits an IOV model with a fixed inter-occasion Omega", {
    .d <- nlmixr2data::theo_sd
    .d$occ <- ifelse(.d$TIME < 5, 1L, 2L)
    .mod <- function() {
      ini({ tka <- log(1.5); tv <- log(31.5); tke <- log(0.08); add.sd <- 0.7
        eta.ka ~ 0.3; eta.ke ~ 0.1
        iov.ka ~ fix(0.05) | occ })
      model({ ka <- exp(tka + eta.ka + iov.ka); v <- exp(tv); ke <- exp(tke + eta.ke)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - ke * center
        cp <- center / v; cp ~ add(add.sd) })
    }
    f <- nlmixr2(.mod, .d, est = "npag",
                 control = npagControl(points = 200L, cycles = 6L))
    expect_s3_class(f, "nlmixr2FitData")
    expect_true(is.finite(as.numeric(f$objf)))
    # the fixed IOV variance is held at 0.05, not estimated from the support points
    expect_equal(unname(f$omega$occ["iov.ka", "iov.ka"]), 0.05, tolerance = 1e-8)
  })

})
