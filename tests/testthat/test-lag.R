nmTest({
  test_that("test lag with warfarin", {
    # Use centralized fit from helper-fits.R
    f <- one.compartment.with.lag.fit.focei

    # objf reflects the current jump-sensitivity solve (~2191)
    expect_true(f$objf < 2300)
  })
})
