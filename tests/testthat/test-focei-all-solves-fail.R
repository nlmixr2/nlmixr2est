nmTest({
  # a fixed-step explicit method on a stiff model (ke0 = 2000/h) fails the
  # ODE solve of every subject; the fit returned the failed-solve penalty
  # (2400 for theo) as its objective instead of an error
  eff <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      tke0 <- fix(7.6)
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      ke0 <- exp(tke0)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      d / dt(ce) <- ke0 * (center / v - ce)
      ce ~ add(add.sd)
    })
  }
  .ctl <- function(method) {
    foceiControl(
      maxOuterIterations = 2,
      print = 0,
      covMethod = "",
      calcTables = FALSE,
      rxControl = rxode2::rxControl(method = method)
    )
  }

  test_that("focei errors when every subject's solve fails at the initial estimates", {
    expect_error(
      suppressMessages(suppressWarnings(nlmixr(eff, nlmixr2data::theo_sd, "focei", .ctl("rk4")))),
      "failed for every subject"
    )
  })

  test_that("the same model fits with a stiff method", {
    .f <- suppressMessages(suppressWarnings(nlmixr(eff, nlmixr2data::theo_sd, "focei", .ctl("liblsoda"))))
    expect_true(is.finite(.f$objf))
    expect_lt(.f$objf, 1000)
  })
})
