nmTest({
  test_that("ini ui works", {
    # Use centralized model from helper-models.R
    f <-  one.compartment

    # fit the model
    fit <- .nlmixr(one.compartment, theo_sd, est = "saem", control = saemControlFast)

    expect_equal(fit$iniUi$iniDf, f$iniDf)
  })
})
