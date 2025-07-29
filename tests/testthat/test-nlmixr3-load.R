nmTest({
  test_that("nlmixr3 fit loads", {

    fit <- readRDS(system.file("testfit_nlmixr3.rds", package = "nlmixr2est"))

    expect_error(rxode2::rxSolve(fit))

    expect_snapshot(
      fit <- suppressMessages(nlmixr2fix(fit))
    )

    expect_error(rxode2::rxSolve(fit), NA)

  })
})
