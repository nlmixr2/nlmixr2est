nmTest({
  test_that("predict", {

    # Use centralized model from helper-models.R
    ## The fit is performed by the function nlmixr/nlmix2 specifying
    ## the model, data and estimate

    fit <- .nlmixr(one.compartment, theo_sd, est = "focei",
                   foceiControl(maxOuterIterations = 0L))

    md <- do.call("predict", c(list(fit, theo_md), fit$control))

    md2 <- .nlmixr(fit, theo_md, "predict")

    expect_equal(as.data.frame(md), as.data.frame(md2), tolerance = 1e-4)

    md <- predict(fit, theo_md)

    expect_equal(as.data.frame(md), as.data.frame(md2), tolerance = 1e-4)

    ipred <- predict(fit, theo_sd, level="individual")
    expect_equal(ipred$ipredSim, fit$IPRED, tolerance = 1e-6)

    # Test explicit population level (default)
    ppred <- predict(fit, theo_sd, level="population")
    expect_true("pred" %fin% names(ppred))

    # Test numeric level=0 (population)
    ppred0 <- predict(fit, theo_sd, level=0)
    expect_true("pred" %fin% names(ppred0))

    # Test numeric level=1 (individual)
    ipred1 <- predict(fit, theo_sd, level=1)
    expect_equal(ipred1$ipredSim, fit$IPRED, tolerance = 1e-6)

    # Test alias level="pred" (population)
    ppredAlias <- predict(fit, theo_sd, level="pred")
    expect_true("pred" %fin% names(ppredAlias))

    # Test alias level="ppred" (population)
    ppredAlias2 <- predict(fit, theo_sd, level="ppred")
    expect_true("pred" %fin% names(ppredAlias2))

    # Test alias level="ipred" (individual)
    ipredAlias <- predict(fit, theo_sd, level="ipred")
    expect_equal(ipredAlias$ipredSim, fit$IPRED, tolerance = 1e-6)

    # Test invalid numeric level throws error
    expect_error(predict(fit, theo_sd, level=2),
                 "level numeric must be 0 \\(population\\) or 1 \\(individual\\)")

    # Test individual predictions with new data (theo_md)
    ipredNewData <- predict(fit, theo_md, level="individual")
    expect_true("ipredSim" %fin% names(ipredNewData))

  })
})
