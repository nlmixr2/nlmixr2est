nlmixrControlTest <- function(control) {
  expect_true(inherits(control$rxControl, "rxControl"))
}

test_that("test control options", {
  expect_error(foceiControl(), NA)
  nlmixrControlTest(foceiControl())
  .ctl <- foceiControl()
  expect_error(do.call(foceiControl, .ctl), NA)
  .ctl2 <- do.call(foceiControl, .ctl)
  expect_equal(.ctl, .ctl2)
})
