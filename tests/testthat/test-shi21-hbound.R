## Regression guard for the shi21 finite-difference step-size bound
## (src/shi21.cpp shi21Central / shi21Forward).
##
## The adaptive Shi (2021) step-size search grows h when it cannot detect
## curvature in the finite differences.  For a FOCEi eta whose finite-difference
## information is (numerically) degenerate the search could grow h without bound,
## perturbing the eta so far that the ODE solve degenerates -- which corrupts the
## shared solver state and collapses the eta finite-difference sensitivity for
## every subsequent inner iteration (the inner problem then leaves the eta stuck
## near 0).  shi21Central/shi21Forward now clamp h to a reasonable region
## (hMin <= h <= hMax) so the probe stays sane.
##
## This exercises the FOCEi FD eta-sensitivity path via a dosing-parameter (lag
## time) eta -- the etas that are always finite-differenced -- and checks the eta
## is genuinely estimated (a collapsed FD sensitivity would leave the lag eta
## unmoved and its EBEs at ~0).

test_that("FOCEi finite-differenced (lag-time) eta is estimated -- shi21 step stays bounded", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  .old <- rxode2::getRxThreads(); on.exit(rxode2::setRxThreads(.old), add = TRUE)
  rxode2::setRxThreads(1L)

  d <- nlmixr2data::theo_sd
  mod <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      tlag <- log(0.15)
      add.sd <- 0.7
      eta.cl ~ 0.1
      eta.lag ~ 0.05
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      tlg <- exp(tlag + eta.lag)      # lag-time IIV -> finite-differenced eta
      lag(depot) <- tlg
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  fit <- suppressWarnings(suppressMessages(
    nlmixr2(mod, d, "focei",
            foceiControl(print = 0L, maxOuterIterations = 20L,
                         maxInnerIterations = 30L, calcTables = FALSE))))

  ## the fit converges to a finite, sane objective
  expect_true(is.finite(fit$objf))
  expect_lt(fit$objf, 260)

  ## the finite-differenced lag eta is genuinely estimated (not collapsed): its
  ## variance is a sane positive value and its EBEs are non-degenerate.
  expect_gt(unname(fit$omega["eta.lag", "eta.lag"]), 0.05)
  expect_gt(diff(range(fit$eta$eta.lag)), 0.5)
})
