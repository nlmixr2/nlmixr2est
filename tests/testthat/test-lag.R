nmTest({
  test_that("test lag with warfarin", {
    # Use centralized fit from helper-fits.R
    f <- one.compartment.with.lag.fit.focei

    expect_true(f$objf < 500)
  })
})
