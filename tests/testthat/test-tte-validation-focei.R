nmTest({
  skip_on_cran()
  skip_if_not(.nlmixr2estRxode2HasTte(), "rxode2 without tte() endpoints")
  skip_if(is.null(.tteValidation()), "tte validation scenarios not installed")
  .ctl <- foceiControl(print = 0)

  test_that("tte validation: exponential matches the closed form", {
    .fit <- .tteValidationCheck("exponential", "bobyqa", bobyqaControl(print = 0), tol = 0.01, llTol = 1e-3)
    expect_equal(
      unname(.tteValidationEst(.fit)["lh"]),
      unname(.tteValidationExactFit("exponential")$est),
      tolerance = 1e-4
    )
  })

  # repeated events determine the frailty well, and the Laplace
  # approximation is close to the exact likelihood
  for (.file in c("rtte", "gapTime", "pkFixed", "tvCov")) {
    test_that(paste("tte validation:", .file, "(focei)"), {
      .fit <- .tteValidationCheck(.file, "focei", .ctl, tol = 0.35, llTol = 2)
      if (.file == "tvCov") {
        .red <- .tteValidationCheck("tvCov-reduced", "focei", .ctl, tol = 0.35, llTol = 2)
        .tteValidationCheckDofv(.fit, .red, "focei", tol = 1)
      }
    })
  }

  # the joint PK/hazard model has no exact likelihood: the truth and any
  # stored reference
  test_that("tte validation: pkJoint (focei)", {
    .fit <- .tteValidationCheck("pkJoint", "focei", .ctl, tol = 0.35)
    .tteValidationExpectTruth(.fit, "pkJoint")
  })

  # a single event per subject barely determines the frailty, and the
  # Laplace approximation is several -2LL units off the exact likelihood
  for (.file in c("weibull", "interval", "lateEntry", "dropout")) {
    test_that(paste("tte validation:", .file, "(focei)"), {
      .tteValidationCheck(.file, "focei", .ctl, tol = 1, llTol = 8)
    })
  }
})
