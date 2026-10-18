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

test_that("a simulation restores every estimation field the reset creates", {
  # .preProcessHooksRun() saves and restores .nlmixr2EstEnvFields around a
  # simulation; a field .nlmixr2globalReset() creates without being listed
  # there would leak from a simulation into the next estimation
  .b <- deparse(body(.nlmixr2globalReset))
  .made <- regmatches(.b, regexpr("(?<=nlmixr2EstEnv\\$)[A-Za-z0-9_.]+(?= <-)", .b, perl = TRUE))
  expect_gt(length(.made), 0L)
  expect_setequal(.made, .nlmixr2EstEnvFields)
})

test_that("every simulation method leaves the hooks' globals as it found them", {
  # .uiApplyIov() clears the IOV rewrite's state and the mu-referencing
  # translations are reset for any method; a simulation inside a fit (its
  # table step's vpcSim() or population predict()) must leave the fit's
  .state <- .preProcessHookState()
  withr::defer(.preProcessHookStateRestore(.state))
  for (.est in .preProcessHooksSimulationEst) {
    .uiIovEnv$ui <- "fit's rewrite"
    .muRefTrans$cur <- "fit's translations"
    .sim <- new.env(parent = emptyenv())
    .sim$ui <- rxode2::assertRxUi(.iovFitMod)
    .sim$data <- .iovFitData()
    suppressMessages(.preProcessHooksRun(.sim, .est))
    expect_identical(.uiIovEnv$ui, "fit's rewrite", label = .est)
    expect_identical(.muRefTrans$cur, "fit's translations", label = .est)
  }
})

# A saem model with a mu-referenced covariate (wt.cl * log(WT / 70)), which
# saem's hooks rewrite into the data and the model
.hookMuCovariate <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    wt.cl <- 0.75
    add.sd <- 0.7
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl + wt.cl * log(WT / 70))
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

nmTest({
  test_that("a saem fit with a mu-referenced covariate asked for NPDE is finished", {
    # npde's vpcSim() in the table step reset the translations the fit's
    # finish undoes
    .fit <- suppressWarnings(suppressMessages(nlmixr2(
      .hookMuCovariate,
      nlmixr2data::theo_sd,
      "saem",
      control = saemControl(print = 0, covMethod = ""),
      table = tableControl(npde = TRUE)
    )))
    expect_true(all(c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD") %in% names(.fit)))
    expect_false(any(grepl("nlmixrMuDerCov", names(.fit))))
    expect_false(any(grepl("nlmixrMuDerCov", names(.fit$dataSav))))
    expect_false(any(grepl("nlmixrMuDerCov", deparse(.fit$ui$fun))))
    expect_identical(.fit$ui$iniDf$name, rxode2::assertRxUi(.hookMuCovariate)$iniDf$name)
  })
})
