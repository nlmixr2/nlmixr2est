nmTest({
  test_that("nlmUnscalePar() is exported (#940)", {
    expect_true("nlmUnscalePar" %in% getNamespaceExports("nlmixr2est"))
  })

  test_that("nlmUnscalePar() refuses when no nlm problem is loaded (#940)", {
    # nlmFree() releases the scale buffers but leaves ntheta at its old value,
    # so without the loaded guard a stale-length vector would read freed
    # memory rather than erroring
    nlmFree()
    expect_error(nlmUnscalePar(c(1, 2)), "not loaded")
    expect_error(nlmScalePar(c(1, 2)), "not loaded")
  })
})
