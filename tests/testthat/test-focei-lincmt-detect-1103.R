nmTest({
  # #1103: the linCmt() guards only checked predDf$linCmt, which is TRUE only when
  # the ENDPOINT is linCmt().  `cp <- linCmt()` slipped through, and under
  # rxode2 >= 5.1.8 its 2nd-order expansion silently drops terms instead of failing,
  # so fast=TRUE gradients and the default analytic covariance were wrong.
  .linLhs <- function() {
    ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
          eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
            cp <- linCmt()
            cp ~ add(add.sd) })
  }

  test_that(".foceiUsesLinCmt sees linCmt() anywhere, not only as the endpoint", {
    .end <- function() {
      ini({ tka <- 0.4; tcl <- 1; tv <- 3.4; eta.ka ~ 0.6; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); linCmt() ~ add(add.sd) })
    }
    .ll <- function() {
      ini({ tka <- 0.4; tcl <- 1; tv <- 3.4; eta.ka ~ 0.6; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt()
              ll(err) ~ -0.5 * log(2 * pi) - log(add.sd) - 0.5 * ((DV - cp) / add.sd)^2 })
    }
    .ode <- function() {
      ini({ tka <- 0.4; tcl <- 1; tv <- 3.4; eta.ka ~ 0.6; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
              d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
              cp <- center / v; cp ~ add(add.sd) })
    }
    expect_true(.foceiUsesLinCmt(rxode2::rxode2(.end)))
    expect_true(.foceiUsesLinCmt(rxode2::rxode2(.linLhs)))
    expect_true(.foceiUsesLinCmt(rxode2::rxode2(.ll)))
    # the endpoint-only flag is what used to miss the last two
    expect_false(isTRUE(any(rxode2::rxode2(.linLhs)$predDf$linCmt)))
    expect_false(.foceiUsesLinCmt(rxode2::rxode2(.ode)))
  })

  test_that("cp <- linCmt(): no analytic covariance, fast=TRUE falls back to the fast=FALSE fit", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    d <- nlmixr2data::theo_sd
    fF <- suppressMessages(suppressWarnings(
      nlmixr2(.linLhs, d, "focei", foceiControl(print = 0L, covMethod = "analytic"))
    ))
    # the analytic covariance declines to the FD route instead of using wrong 2nd derivatives
    expect_equal(fF$covMethod, "r,s (full)")
    fT <- suppressMessages(suppressWarnings(
      nlmixr2(.linLhs, d, "focei", foceiControl(print = 0L, covMethod = "", fast = TRUE))
    ))
    expect_false(fT$control$fast)
    # the downgrade re-defaults the outer optimizer too (lbfgsb3c stalled at 133.57)
    expect_equal(fT$control$outerOptTxt, "bobyqa")
    expect_equal(fT$objf, fF$objf, tolerance = 1e-4)

    # ground truth: the ODE twin's covariance, same estimator, at the same estimates
    .ode <- function() {
      ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
            eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              d/dt(depot) <- -ka * depot
              d/dt(center) <- ka * depot - cl / v * center
              cp <- center / v
              cp ~ add(add.sd) })
    }
    ui <- rxode2::rxode2(.ode)
    iniDf <- ui$iniDf
    iniDf$est <- fF$finalUi$iniDf$est[match(iniDf$name, fF$finalUi$iniDf$name)]
    ui$iniDf <- iniDf
    tight <- rxode2::rxControl(atol = 1e-10, rtol = 1e-10, atolSens = 1e-10, rtolSens = 1e-10)
    fO <- suppressMessages(suppressWarnings(nlmixr2(
      ui,
      d,
      "focei",
      foceiControl(
        print = 0L,
        covMethod = "r,s",
        maxOuterIterations = 0L,
        rxControl = tight
      )
    )))
    expect_equal(fO$covMethod, "r,s (full)")
    # the old linCmt() analytic SE was up to 37% off; the fallback matches the twin
    expect_equal(sqrt(diag(fF$cov)), sqrt(diag(fO$cov)), tolerance = 0.02)
  })

  test_that("full conditional Hessian refuses linCmt() with a clear error", {
    skip_if_not_installed("nlmixr2data")
    ctl <- list(maxOuterIterations = 0L, covMethod = "", calcTables = FALSE, print = 0L)
    expect_error(
      suppressMessages(nlmixr2(.linLhs, nlmixr2data::theo_sd, est = "flaplace", control = ctl)),
      "linCmt"
    )
  })
})
