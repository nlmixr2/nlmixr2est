nmTest({
  skip_on_cran()
  skip_if_not(.nlmixr2estRxode2HasTte(), "rxode2 without tte() endpoints")
  skip_if(is.null(.tteValidation()), "tte validation scenarios not installed")
  .ctl <- saemControl(print = 0)

  # saem is compared with Monolix's SAEM.  Its objective is a 3-node
  # quadrature, so only its estimates are compared.  It is not run where the solve has extra inputs (a data
  # covariate or the CMT of a second endpoint: pkFixed, tvCov, dropout), where
  # it currently diverges, nor on the slowest scenarios.
  test_that("tte validation: rtte (saem)", {
    .tteValidationCheck("rtte", "saem", .ctl, tol = 0.5, ref = "monolix", refTol = 1.5)
  })
  test_that("tte validation: lateEntry (saem)", {
    .tteValidationCheck("lateEntry", "saem", .ctl, tol = 1, ref = "monolix", refTol = 2)
  })
})
