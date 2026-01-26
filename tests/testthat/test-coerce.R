nmTest({
  # Use centralized fit from helper-fits.R
  fit <- one.compartment.fit.saem

  test_that("as.rxUi works for estimated models", {
    expect_s3_class(as.rxUi(fit), "rxUi")
  })
})
