nlmixrControlTest <- function(control) {
  expect_true(inherits(control$rxControl, "rxControl"))
  expect_true(checkmate::testNumeric(control$ci, lower=0, upper=1, any.missing=FALSE, len=1))
  expect_true(checkmate::testIntegerish(control$sigdigTable, lower=1, any.missing=FALSE, len=1))
}

test_that("test foceiControl option sanity", {
  expect_error(foceiControl(), NA)
  nlmixrControlTest(foceiControl())
  .ctl <- foceiControl()
  expect_error(do.call(foceiControl, .ctl), NA)
  .ctl2 <- do.call(foceiControl, .ctl)
  expect_equal(.ctl, .ctl2)

  # ResetEtaP
  .ctl <- foceiControl(resetEtaP=0.5)
  .ctl2 <- do.call(foceiControl, .ctl)
  expect_equal(.ctl, .ctl2)
  # resetThetaP
  .ctl <- foceiControl(resetThetaP=0.5)
  .ctl2 <- do.call(foceiControl, .ctl)
  expect_equal(.ctl, .ctl2)
  # resetThetaFinalP
  .ctl <- foceiControl(resetThetaFinalP=0.5)
  .ctl2 <- do.call(foceiControl, .ctl)
  expect_equal(.ctl, .ctl2)

})
