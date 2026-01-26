nmTest({
  # Use centralized model from helper-models.R
  fit <- .nlmixr(one.compartment, theo_sd, est = "saem", control = saemControlFast)

  test_that("as.rxUi works for estimated models", {
    expect_s3_class(as.rxUi(fit), "rxUi")
  })
})
