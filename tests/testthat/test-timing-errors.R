# A fit whose session information cannot be read fails after its output
# environment is made the timing environment, and before it has a time table
.timingSessionInfoFails <- function() stop("session information failed on purpose")

.timingOneCmt <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
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

nmTest({
  test_that("the overall timing skips a timing environment with no time table", {
    withr::defer(.popNlmixr2Timing())
    nlmixr2global$currentTimingEnvironment <- new.env(parent = emptyenv())
    nlmixr2global$nlmixr2Time <- proc.time()
    expect_no_error(.finalizeOverallTiming())
    expect_null(nlmixr2global$currentTimingEnvironment)
    expect_null(nlmixr2global$nlmixr2Time)
  })

  test_that("a timed stage keeps its code's error, and its time, without a time table", {
    withr::defer(.popNlmixr2Timing())
    .env <- new.env(parent = emptyenv())
    nlmixr2global$currentTimingEnvironment <- .env
    expect_identical(
      tryCatch(nlmixrWithTiming("stage", stop("the stage failed")), error = conditionMessage),
      "the stage failed"
    )
    # the stage goes to the extra table, as with no timing environment
    expect_false(exists("time", envir = .env, inherits = FALSE))
    expect_named(nlmixr2global$extraTimingTable, "stage")
    expect_null(nlmixr2global$timingStack)
  })

  test_that("a fit that fails after its output is started reports its own error", {
    skip_on_cran()
    local_mocked_bindings(.sessionInfo = .timingSessionInfoFails)
    .ctl <- list(
      saem = saemControl(print = 0, nBurn = 5, nEm = 5, covMethod = "", calcTables = FALSE),
      focei = foceiControl(print = 0, maxOuterIterations = 0, covMethod = "", calcTables = FALSE)
    )
    for (.est in names(.ctl)) {
      .e <- tryCatch(
        suppressMessages(nlmixr2(.timingOneCmt, nlmixr2data::theo_sd, .est, control = .ctl[[.est]])),
        error = conditionMessage
      )
      expect_identical(.e, "session information failed on purpose", label = .est)
    }
    # the timing state is reset for the next fit
    expect_null(nlmixr2global$currentTimingEnvironment)
  })
})
