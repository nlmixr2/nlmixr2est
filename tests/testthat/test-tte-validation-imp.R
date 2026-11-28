nmTest({
  skip_on_cran()
  skip_if_not(.nlmixr2estRxode2HasTte(), "rxode2 without tte() endpoints")
  skip_if(is.null(.tteValidation()), "tte validation scenarios not installed")
  .ctl <- impControl(print = 0)

  # compared with Monolix.  imp's reported objective is the Laplace one at
  # its estimates, so only its estimates are compared; the slowest scenarios are left to focei and agq
  for (.file in c("weibull", "rtte", "pkFixed", "tvCov", "dropout")) {
    test_that(paste("tte validation:", .file, "(imp)"), {
      .fit <- .tteValidationCheck(.file, "imp", .ctl, tol = 0.5, ref = "monolix", refTol = 1.5)
      if (.file == "tvCov") {
        # the time-varying covariate stays out of the baseline regression
        expect_lt(.tteValidationEst(.fit)[["bx"]], -0.3)
      }
    })
  }
})
