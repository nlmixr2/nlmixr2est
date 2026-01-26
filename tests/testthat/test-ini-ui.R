nmTest({
  test_that("ini ui works", {
    # Use centralized model from helper-models.R
    f <-  one.compartment

    # Use centralized fit from helper-fits.R
    fit <- one.compartment.fit.saem

    expect_equal(fit$iniUi$iniDf, f$iniDf)
  })
})
