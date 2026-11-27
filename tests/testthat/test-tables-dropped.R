nmTest({
  # #1354: a failed tables step used to leave only a $runInfo line
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
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
  .ctl <- foceiControl(print = 0L, maxOuterIterations = 0L, covMethod = "")

  test_that("a dropped focei tables step warns and flags the fit", {
    skip_on_cran()
    local_mocked_bindings(addTable = function(...) stop("table boom"))
    .w <- character(0)
    fit <- withCallingHandlers(
      suppressMessages(nlmixr(.mod, theo_sd, "focei", .ctl)),
      warning = function(w) {
        .w <<- c(.w, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_true(any(grepl("tables step failed", .w, fixed = TRUE)))
    expect_true(fit$tablesDropped)
    expect_false(inherits(fit, "data.frame"))
    expect_true(any(grepl("returning without table step", fit$runInfo, fixed = TRUE)))
  })

  test_that("a fit with its tables is not flagged", {
    skip_on_cran()
    fit <- suppressMessages(.nlmixr(.mod, theo_sd, "focei", .ctl))
    expect_false(fit$tablesDropped)
    expect_false(any(grepl("tables step failed", fit$runInfo, fixed = TRUE)))
    # tables turned off on purpose are not dropped
    .ctl2 <- .ctl
    .ctl2$calcTables <- FALSE
    fit2 <- suppressMessages(.nlmixr(.mod, theo_sd, "focei", .ctl2))
    expect_false(fit2$tablesDropped)
  })

  test_that("adding the tables afterwards clears the flag", {
    .env <- new.env(parent = emptyenv())
    .tablesDroppedSet(.env, TRUE)
    expect_true(.tablesDroppedGet(.env))
    .df <- data.frame(a = 1)
    .cls <- c("nlmixr2FitData", "data.frame")
    attr(.cls, ".foceiEnv") <- .env
    class(.df) <- .cls
    expect_true(.tablesDroppedGet(.df))
    .tablesDroppedSet(.df, FALSE)
    expect_false(.tablesDroppedGet(.env))
    expect_false(.tablesDroppedGet(new.env(parent = emptyenv())))
    expect_false(.tablesDroppedGet("fit"))
  })

  test_that("a degraded saem output flags the fit", {
    .fitEnv <- new.env(parent = emptyenv())
    .n <- 0L
    local_mocked_bindings(nlmixr2CreateOutputFromUi = function(ui, data, control, table, env, est) {
      .n <<- .n + 1L
      if (.n == 1L) {
        stop("build failure")
      }
      .fitEnv
    })
    .e <- new.env(parent = emptyenv())
    .e$control <- list(calcTables = TRUE)
    .e$table <- list(cwres = TRUE, npde = TRUE)
    expect_warning(.saemCreateOutput(.e), "fit degraded, no tables")
    expect_true(.tablesDroppedGet(.fitEnv))
  })
})
