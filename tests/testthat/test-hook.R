# The pre-processing hooks (R/hook.R).

.hookFixedTheta <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- fix(3.45)
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    linCmt() ~ add(add.sd)
  })
}

test_that("a simulation's pre-processing leaves the estimation state as it found it", {
  .old <- mget(.nlmixr2EstEnvFields, envir = nlmixr2global$nlmixr2EstEnv, ifnotfound = list(NULL))
  withr::defer(list2env(.old, envir = nlmixr2global$nlmixr2EstEnv))
  nlmixr2global$nlmixr2EstEnv$uiUnfix <- NULL
  # augPred(), vpcSim() and an nlmixr2() refit's simulation information run the
  # hooks with no control, so the fixed theta is substituted into the model ...
  .sim <- new.env(parent = emptyenv())
  .sim$ui <- rxode2::assertRxUi(.hookFixedTheta)
  .sim$data <- nlmixr2data::theo_sd
  suppressMessages(.preProcessHooksRun(.sim, "rxSolve"))
  expect_false("tv" %in% .sim$ui$iniDf$name)
  # ... without recording the unsubstituted model for an estimation to restore
  expect_null(nlmixr2global$nlmixr2EstEnv$uiUnfix)
  # an estimation that substitutes it does record it
  .est <- new.env(parent = emptyenv())
  .est$ui <- rxode2::assertRxUi(.hookFixedTheta)
  .est$data <- nlmixr2data::theo_sd
  .est$control <- foceiControl(literalFix = TRUE)
  suppressMessages(.preProcessHooksRun(.est, "focei"))
  expect_true("tv" %in% nlmixr2global$nlmixr2EstEnv$uiUnfix$iniDf$name)
  # and a simulation run afterwards keeps it
  suppressMessages(.preProcessHooksRun(.sim, "rxSolve"))
  expect_true("tv" %in% nlmixr2global$nlmixr2EstEnv$uiUnfix$iniDf$name)
})

nmTest({
  test_that("a fit that keeps its fixed theta in the model refits, after a simulation too", {
    .ctl <- foceiControl(print = 0L, literalFix = FALSE, covMethod = "")
    .fit <- suppressMessages(nlmixr2(.hookFixedTheta, nlmixr2data::theo_sd, "focei", control = .ctl))
    expect_identical(unname(.fit$parFixedDf["tv", "Estimate"]), 3.45)
    # nlmixr2(<fit>) builds the fit's simulation information before it refits
    .refit <- suppressMessages(nlmixr2(.fit, nlmixr2data::theo_sd, est = "focei", control = .ctl))
    expect_s3_class(.refit, "nlmixr2FitData")
    expect_identical(unname(.refit$parFixedDf["tv", "Estimate"]), 3.45)
    expect_true(is.na(.refit$parFixedDf["tv", "SE"]))
    # augPred() simulates; a covariance installed afterwards rebuilds the table
    expect_s3_class(suppressMessages(augPred(.fit)), "data.frame")
    suppressWarnings(suppressMessages(setCov(.fit, "r")))
    expect_identical(.fit$covMethod, "r")
    expect_identical(unname(.fit$parFixedDf["tv", "Estimate"]), 3.45)
  })
})
