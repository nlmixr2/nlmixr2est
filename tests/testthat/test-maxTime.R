test_that("the time limit requests one stop and is reset at the end", {
  on.exit({
    .nmStopFileSet("")
    .nmStopTimeSet(NA_real_)
  })
  expect_false(.nmStopActive())
  .nmStopTimeSet(3600)
  expect_true(.nmStopActive())
  expect_false(.nmStopFilePoll())
  .nmStopTimeSet(0)
  expect_true(.nmStopFilePoll())
  expect_true(.nmStopFileRequested())
  expect_true(.nmStopTimedOut())
  expect_equal(.nmStopMessage(), "stopped at the time limit")
  # a request is consumed by the loop it ends; the limit does not fire again
  .nmStopFileSet("")
  expect_false(.nmStopFilePoll())
  expect_true(.nmStopActive())
  expect_warning(expect_true(.nmStopFileEnd()), "stopped at the time limit")
  expect_false(.nmStopActive())
  expect_false(.nmStopTimedOut())
  expect_equal(.nmStopMessage(), "stopped early by user request")
})

test_that("a stop file found with the time up is removed and reported", {
  .sf <- tempfile("stop")
  on.exit({
    .nmStopFileSet("")
    .nmStopTimeSet(NA_real_)
    unlink(.sf)
  })
  expect_true(.nmStopFileStart(list(stopFile = .sf, maxTime = 0)))
  file.create(.sf)
  expect_true(.nmStopFilePoll())
  expect_false(file.exists(.sf))
  expect_false(.nmStopTimedOut())
  expect_warning(expect_true(.nmStopFileEnd()), "stopped early by user request")
})

test_that("a stop file after a time-out is reported as a user request", {
  .sf <- tempfile("stop")
  on.exit({
    .nmStopFileSet("")
    .nmStopTimeSet(NA_real_)
    unlink(.sf)
  })
  expect_true(.nmStopFileStart(list(stopFile = .sf, maxTime = 0)))
  expect_true(.nmStopFilePoll())
  expect_equal(.nmStopMessage(), "stopped at the time limit")
  # the time-out was handled; the file is still polled
  .nmStopFileSet(.sf)
  file.create(.sf)
  expect_true(.nmStopFilePoll())
  expect_equal(.nmStopMessage(), "stopped early by user request")
  .w <- character(0)
  withCallingHandlers(
    .nmStopFileEnd(),
    warning = function(w) {
      .w <<- c(.w, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_setequal(.w, c("stopped at the time limit", "stopped early by user request"))
})

test_that("a nested estimation keeps the outer time limit", {
  on.exit({
    .nmStopFileSet("")
    .nmStopTimeSet(NA_real_)
  })
  expect_false(.nmStopFileStart(list(maxTime = NULL)))
  expect_true(.nmStopFileStart(list(maxTime = 0)))
  expect_false(.nmStopFileStart(list(maxTime = 3600)))
  expect_true(.nmStopFilePoll())
  expect_warning(expect_true(.nmStopFileEnd()), "stopped at the time limit")
  expect_false(.nmStopActive())
})

test_that("foceiControl(maxTime=) is validated", {
  expect_null(foceiControl()$maxTime)
  expect_equal(foceiControl(maxTime = 10)$maxTime, 10)
  expect_error(foceiControl(maxTime = -1))
  expect_error(foceiControl(maxTime = Inf))
  expect_error(foceiControl(maxTime = "1"))
})

skip_on_cran()

.maxTimeMod <- function() {
  ini({
    tka <- 0.45
    tcl <- log(c(0, 2.7, 100))
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

test_that("maxTime ends a FOCEi fit at the time limit", {
  .fit <- suppressMessages(nlmixr(
    .maxTimeMod,
    nlmixr2data::theo_sd,
    "focei",
    control = foceiControl(print = 0, maxTime = 0, covMethod = "")
  ))
  expect_equal(.fit$message, "stopped at the time limit")
  expect_true("stopped at the time limit" %in% .fit$runInfo)
  expect_false("stopped early by user request" %in% .fit$runInfo)
  # stopped at the first evaluation: the initial estimates
  expect_equal(.fit$theta[["tka"]], 0.45, tolerance = 1e-6)
  .ph <- .fit$parHistData
  expect_length(.ph$objf[.ph$type == "Unscaled"], 1L)
  expect_false(.nmStopActive())

  .fit2 <- suppressMessages(nlmixr(
    .maxTimeMod,
    nlmixr2data::theo_sd,
    "focei",
    control = foceiControl(print = 0, maxTime = NULL, covMethod = "")
  ))
  expect_false(.fit2$message == "stopped at the time limit")
  expect_false("stopped at the time limit" %in% .fit2$runInfo)
  .ph2 <- .fit2$parHistData
  expect_gt(length(.ph2$objf[.ph2$type == "Unscaled"]), 1L)
})

test_that("a generous maxTime does not stop a FOCEi fit", {
  .fit <- suppressMessages(nlmixr(
    .maxTimeMod,
    nlmixr2data::theo_sd,
    "focei",
    control = foceiControl(print = 0, maxTime = 1e6, covMethod = "")
  ))
  expect_false("stopped at the time limit" %in% .fit$runInfo)
})

test_that("maxTime ends a SAEM fit at the current iteration", {
  expect_null(saemControl()$maxTime)
  expect_equal(saemControl(maxTime = 10)$maxTime, 10)
  expect_error(saemControl(maxTime = -1))
  .fit <- suppressMessages(nlmixr(
    .maxTimeMod,
    nlmixr2data::theo_sd,
    "saem",
    control = saemControl(print = 0, maxTime = 0, nBurn = 20, nEm = 30)
  ))
  expect_true("stopped at the time limit" %in% .fit$runInfo)
  expect_equal(nrow(.fit$parHist), 1L)
  expect_true(is.finite(.fit$objf))
  expect_false(.nmStopActive())

  .fit2 <- suppressMessages(nlmixr(
    .maxTimeMod,
    nlmixr2data::theo_sd,
    "saem",
    control = saemControl(print = 0, maxTime = 1e6, nBurn = 20, nEm = 30)
  ))
  expect_false("stopped at the time limit" %in% .fit2$runInfo)
  expect_equal(nrow(.fit2$parHist), 50L)
})

.maxTimeModNlm <- function() {
  ini({
    tka <- 0.45
    tcl <- log(c(0, 2.7, 100))
    tv <- 3.45
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl)
    v <- exp(tv)
    linCmt() ~ add(add.sd)
  })
}

for (.est in c("nlm", "nlminb", "bobyqa", "uobyqa", "newuoa", "lbfgsb3c", "n1qn1", "optim", "trust")) {
  test_that(paste0("maxTime ends an est='", .est, "' fit early"), {
    .ctl <- paste0(.est, "Control")
    expect_null(do.call(.ctl, list())$maxTime)
    expect_equal(do.call(.ctl, list(maxTime = 10))$maxTime, 10)
    expect_error(do.call(.ctl, list(maxTime = -1)))
    .fit <- suppressMessages(nlmixr(
      .maxTimeModNlm,
      nlmixr2data::theo_sd,
      .est,
      control = list(print = 0, maxTime = 0)
    ))
    expect_equal(.fit$message, "stopped at the time limit")
    expect_true("stopped at the time limit" %in% .fit$runInfo)
    expect_equal(.fit$theta[["tka"]], 0.45, tolerance = 1e-6)
    expect_true(is.finite(.fit$objf))
    expect_false(.nmStopActive())

    .fit2 <- suppressMessages(nlmixr(
      .maxTimeModNlm,
      nlmixr2data::theo_sd,
      .est,
      control = list(print = 0, maxTime = 1e6)
    ))
    expect_false("stopped at the time limit" %in% .fit2$runInfo)
    expect_lt(.fit2$objf, .fit$objf)
  })
}
