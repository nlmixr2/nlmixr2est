nmTest({
  test_that("predict", {
    # Use centralized fit from helper-fits.R
    ## The fit was performed with focei, maxOuterIterations = 0L, and print = 0

    fit <- one.compartment.fit.focei

    # given no solving options both entry points inherit the fit's own rxControl
    # (rtol=1e-3); pass it explicitly here so the pairing is visible
    md <- suppressMessages(do.call("predict", c(list(fit, theo_md), fit$control$rxControl)))

    md2 <- .nlmixr(fit, theo_md, "predict")

    expect_equal(as.data.frame(md), as.data.frame(md2), tolerance = 1e-4)

    # with no solving options both entry points use the fit's own rxControl
    expect_equal(
      as.data.frame(suppressMessages(predict(fit, theo_md))),
      as.data.frame(.nlmixr(fit, theo_md, "predict")),
      tolerance = 1e-4
    )

    # so comparing against rxControl() defaults needs them passed on both sides
    md <- suppressMessages(
      do.call("predict", c(list(fit, theo_md), rxode2::rxControl()))
    )

    md2 <- .nlmixr(fit, theo_md, "predict", control = rxode2::rxControl())

    expect_equal(as.data.frame(md), as.data.frame(md2), tolerance = 1e-4)

    ipred <- suppressMessages(predict(fit, theo_sd, level = "individual"))
    expect_equal(ipred$ipredSim, fit$IPRED, tolerance = 1e-4)

    # Test explicit population level (default)
    ppred <- suppressMessages(predict(fit, theo_sd, level = "population"))
    expect_true("pred" %in% names(ppred))

    # Test numeric level=0 (population)
    ppred0 <- suppressMessages(predict(fit, theo_sd, level = 0))
    expect_true("pred" %in% names(ppred0))

    # Test numeric level=1 (individual)
    ipred1 <- suppressMessages(predict(fit, theo_sd, level = 1))
    expect_equal(ipred1$ipredSim, fit$IPRED, tolerance = 1e-4)

    # Test alias level="pred" (population)
    ppredAlias <- suppressMessages(predict(fit, theo_sd, level = "pred"))
    expect_true("pred" %in% names(ppredAlias))

    # Test alias level="ppred" (population)
    ppredAlias2 <- suppressMessages(predict(fit, theo_sd, level = "ppred"))
    expect_true("pred" %in% names(ppredAlias2))

    # Test alias level="ipred" (individual)
    ipredAlias <- suppressMessages(predict(fit, theo_sd, level = "ipred"))
    expect_equal(ipredAlias$ipredSim, fit$IPRED, tolerance = 1e-4)

    # Test invalid numeric level throws error
    expect_error(
      suppressMessages(predict(fit, theo_sd, level = 2)),
      "level numeric must be 0 \\(population\\) or 1 \\(individual\\)"
    )

    # Test individual predictions with new data (theo_md)
    ipredNewData <- suppressMessages(predict(fit, theo_md, level = "individual"))
    expect_true("ipredSim" %in% names(ipredNewData))
  })
})
