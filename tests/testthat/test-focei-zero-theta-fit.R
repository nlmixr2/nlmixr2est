# A theta that starts at exactly 0 is estimated like one started elsewhere (#1228).
nmTest({
  .d1228 <- function() {
    d <- nlmixr2data::theo_sd
    d$COV <- d$ID %% 2
    obs <- d$EVID == 0 & d$COV == 1
    d$DV[obs] <- d$DV[obs] * 1.5
    d
  }
  .m1228 <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; e_cov_v <- 0
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v + e_cov_v * COV)
      linCmt() ~ add(add.sd)
    })
  }
  .ctl <- function(...) foceiControl(print = 0, calcTables = FALSE, covMethod = "", ...)

  test_that("a covariate effect started at 0 is estimated (#1228)", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    d <- .d1228()
    fit0 <- suppressMessages(nlmixr2(.m1228, d, est = "focei", control = .ctl()))
    fit01 <- suppressMessages(nlmixr2(.m1228 |> ini(e_cov_v = 0.1), d, est = "focei", control = .ctl()))
    # truth is log(1/1.5) = -0.405; frozen it stayed at ~0 with objf ~12 higher
    expect_equal(fit0$theta[["e_cov_v"]], fit01$theta[["e_cov_v"]], tolerance = 0.05)
    expect_lt(fit0$objf, fit01$objf + 0.5)
    # the one-step start moved it most of the way before the optimizer started
    expect_true(fit0$env$zeroThetaStart < -0.3)
    # the typical-scale fix alone also estimates it
    fitS <- suppressMessages(nlmixr2(.m1228, d, est = "focei", control = .ctl(zeroThetaStep = 0)))
    expect_equal(fitS$theta[["e_cov_v"]], fit01$theta[["e_cov_v"]], tolerance = 0.05)
  })
})
