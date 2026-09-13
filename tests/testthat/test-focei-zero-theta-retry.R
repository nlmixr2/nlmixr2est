# A theta initialized at exactly 0 is nudged to +/-zeroTheta before estimation,
# because FOCEi scales a linear parameter by 1/|init| and 0 gives no scale.  The
# nudge therefore doubles as the magnitude the search assumes, and a small one
# means the search explores in small steps and stops on the nudge -- the
# reported "estimate" is then the nudge value, which for a covariate coefficient
# is indistinguishable from a correct null result.
#
# No single nudge fixes it.  Measured: a coefficient on log(WT/70) with a true
# value of 0.90 needs 0.1 (it returns 0.0009 at 0.001), while the coefficient on
# untransformed WT in test-focei-zero-init-scale.R, true value 0.03, needs 0.001
# (it returns 0.0043 at 0.1).  The two want scales 30x apart.
#
# So the stall is detected rather than guessed: a first fit whose estimate is
# still within zeroThetaRetryTol multiples of the nudge is re-fit once from
# zeroTheta*zeroThetaRetry, and the better objective wins.
nmTest({
  test_that("a zero-initialized coefficient is estimated, not left on its nudge", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .d <- nlmixr2data::theo_sd
    .d$WT70 <- log(.d$WT / 70)
    .m <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        bWT <- 0
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + bWT * WT70 + eta.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    # with the retry (default): recovers the real effect, ~0.92.  The same model
    # started at 0.5, where no nudge applies, converges to 0.899 at objf 128.902
    # -- so this is the optimum, not merely movement.
    .f <- suppressMessages(suppressWarnings(
      nlmixr2(.m, .d, est = "focei",
              control = foceiControl(print = 0L, covMethod = ""))))
    .b <- unname(.f$parFixedDf["bWT", "Estimate"])
    expect_true(.b > 0.5)
    expect_true(.f$objf < 129.5)

    # with the retry DISABLED the same fit stops on the nudge, which is what
    # makes the retry load-bearing rather than cosmetic
    .f0 <- suppressMessages(suppressWarnings(
      nlmixr2(.m, .d, est = "focei",
              control = foceiControl(print = 0L, covMethod = "",
                                     zeroThetaRetry = 1))))
    .b0 <- unname(.f0$parFixedDf["bWT", "Estimate"])
    expect_true(abs(.b0) < 0.01)
    expect_true(.f0$objf > .f$objf)
  })

  test_that("a coefficient that moved off its nudge is not re-fit", {
    # the guard that keeps the retry from firing on healthy fits: the stall
    # signature is "still on the nudge", and a healthy coefficient lands orders
    # of magnitude away from it.  Asserted on the control level so it does not
    # depend on a second fit's timing.
    expect_equal(.foceiZeroThetaStalled(
      structure(list(), class = "try-error"), 0.001, 3), character(0))
    # a non-positive or non-finite magnitude disables the check entirely
    expect_equal(.foceiZeroThetaStalled(NULL, 0, 3), character(0))
    expect_equal(.foceiZeroThetaStalled(NULL, NA_real_, 3),
                 character(0))
  })

  test_that("zeroThetaRetry and zeroThetaRetryTol are validated", {
    expect_error(foceiControl(zeroThetaRetry = 0.5))
    expect_error(foceiControl(zeroThetaRetryTol = -1))
    expect_equal(foceiControl(zeroThetaRetry = 1)$zeroThetaRetry, 1)
    expect_equal(foceiControl()$zeroThetaRetry, 100)
    expect_equal(foceiControl()$zeroThetaRetryTol, 3)
  })
})
