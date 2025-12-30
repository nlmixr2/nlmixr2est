nmTest({
  test_that("nlmixr3 fit loads", {
    skip_if_not_installed("qs") # nlmixr3 uses qs for serialization

    fit <- readRDS(system.file("testfit_nlmixr3.rds", package = "nlmixr2est"))

    expect_error(rxode2::rxSolve(fit))

    fit <- expect_error(
      suppressWarnings(suppressMessages(nlmixr2fix(fit))), NA)

    expect_error(rxode2::rxSolve(fit), NA)

  })
})
