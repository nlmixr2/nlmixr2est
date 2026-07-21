# foceiControl(shi21hMax=, shi21hMin=) expose the adaptive shi21 finite-
# difference step bounds and the C++ FD path must actually use them.
nmTest({
  test_that("shi21hMax/shi21hMin round-trip and are validated", {
    ctl <- foceiControl()
    expect_equal(ctl$shi21hMax, 2.0)
    expect_equal(ctl$shi21hMin, 1e-4)
    ctl2 <- foceiControl(shi21hMax = 5, shi21hMin = 1e-5)
    expect_equal(ctl2$shi21hMax, 5)
    expect_equal(ctl2$shi21hMin, 1e-5)
    # hMax must exceed hMin, and both must be finite and non-negative
    expect_error(foceiControl(shi21hMax = 1e-6, shi21hMin = 1e-3), "greater than")
    expect_error(foceiControl(shi21hMax = -1))
    expect_error(foceiControl(shi21hMax = Inf))
  })

  test_that("the shi21 FD step bounds reach the C++ outer gradient", {
    skip_on_cran()
    ode.cmt <- function() {
      ini({
        tka <- 0.3; tcl <- 0.9; tv <- 3.3; add.sd <- 0.7
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.sd)
      })
    }
    # route the outer theta gradient through shi21 (shi21maxOuter != 0) and use a
    # gradient-based outer optimizer so the bounds are actually exercised
    run <- function(...) {
      f <- suppressWarnings(suppressMessages(
        nlmixr2(ode.cmt, theo_sd, "focei",
                foceiControl(print = 0L, covMethod = "", calcTables = FALSE,
                             outerOpt = "nlminb", shi21maxOuter = 20L, ...))))
      unname(f$parFixedDf$Estimate)
    }
    wide <- run(shi21hMax = 2.0, shi21hMin = 1e-4)
    # collapse the step into the roundoff/solver-noise regime -> a materially
    # different (degraded) outer gradient, hence a different converged estimate
    noisy <- run(shi21hMax = 1e-7, shi21hMin = 1e-8)
    expect_false(isTRUE(all.equal(wide, noisy, tolerance = 1e-4)))
  })
})
