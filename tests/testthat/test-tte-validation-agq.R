nmTest({
  skip_on_cran()
  skip_if_not(.nlmixr2estRxode2HasTte(), "rxode2 without tte() endpoints")
  skip_if(is.null(.tteValidation()), "tte validation scenarios not installed")
  # the default 3 nodes are several -2LL units off the exact likelihood for
  # a single-event frailty; 9 nodes reproduce it
  .ctl <- agqControl(print = 0, nAGQ = 9)

  for (.file in c("weibull", "rtte", "gapTime", "pkFixed", "tvCov", "interval", "lateEntry", "dropout")) {
    test_that(paste("tte validation:", .file, "(agq)"), {
      .fit <- .tteValidationCheck(.file, "agq", .ctl, tol = 0.25, llTol = 0.5)
      if (.file == "tvCov") {
        .red <- .tteValidationCheck("tvCov-reduced", "agq", .ctl, tol = 0.25, llTol = 0.5)
        .tteValidationCheckDofv(.fit, .red, "agq", tol = 0.5)
      }
    })
  }
})
