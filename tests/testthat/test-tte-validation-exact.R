nmTest({
  skip_on_cran()
  skip_if_not(.nlmixr2estRxode2HasTte(), "rxode2 without tte() endpoints")
  skip_if(is.null(.tteValidation()), "tte validation scenarios not installed")

  test_that("the exact likelihood of the simulated data recovers the truth", {
    # checks rxode2's tte() simulation against a likelihood independent of it
    .v <- .tteValidation()
    for (.name in names(.v$env$tteValidationScenarios())) {
      .ex <- .tteValidationExactFit(.name)
      if (is.null(.ex)) {
        next
      }
      .ini <- suppressMessages(rxode2::rxode2(.v$env$tteValidationScenarios(.name)[[1]]$model))$iniDf
      .truth <- stats::setNames(.ini$est, .ini$name)
      for (.p in names(.ex$est)) {
        expect_lt(abs(.ex$est[[.p]] - .truth[[.p]]), 3.5 * .ex$se[[.p]], label = paste(.name, .p, "vs truth / se"))
      }
    }
  })

  test_that("the exact exponential fit is the closed form", {
    .ex <- .tteValidationExactFit("exponential")
    .d <- utils::read.csv(file.path(.tteValidation()$dir, "data", "exponential.csv"))
    .last <- .d[!duplicated(.d$ID, fromLast = TRUE), ]
    expect_equal(unname(.ex$est), log(sum(.last$DV) / sum(.last$TIME)))
  })
})
