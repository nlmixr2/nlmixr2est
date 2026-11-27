# saemControl(autoStop=TRUE): noise-aware windowed early stop of the EM phase (#1411)
nmTest({
  test_that("a flat trajectory stops exactly when the window fills and the count is met", {
    .r <- saemAutoStopTest_(matrix(c(1, -2, 0.3), 200, 3, byrow = TRUE), 50L, 4L, 1e-3, 1e-4)
    # stationary windows end at rows 50, 51, 52, 53
    expect_equal(.r$iter, 53L)
    expect_equal(.r$drift, c(0, 0, 0))
    .r <- saemAutoStopTest_(matrix(1, 200, 1), 20L, 1L, 1e-3, 1e-4)
    expect_equal(.r$iter, 20L)
  })

  test_that("the stop never fires before the window fills", {
    .r <- saemAutoStopTest_(matrix(1, 49, 2), 50L, 1L, 1e-3, 1e-4)
    expect_true(is.na(.r$iter))
    expect_true(all(is.na(.r$drift)))
  })

  test_that("a trending trajectory does not stop", {
    .t <- seq(0, 1, length.out = 400)
    .r <- saemAutoStopTest_(cbind(.t, 1), 50L, 4L, 1e-3, 1e-4)
    expect_true(is.na(.r$iter))
    # a trend buried in noise is still a trend once the window spans it
    .noisy <- cbind(.t + rep(c(-0.002, 0.002), 200), 2)
    expect_true(is.na(saemAutoStopTest_(.noisy, 50L, 4L, 1e-3, 1e-4)$iter))
  })

  test_that("a noisy stationary trajectory stops; the Monte-Carlo error sets the tolerance", {
    # deterministic pseudo-noise around a constant; rtol alone would never pass
    .n <- 300
    .e <- sin(seq_len(.n) * 2.3) * 0.2 + cos(seq_len(.n) * 0.7) * 0.1
    .r <- saemAutoStopTest_(cbind(5 + .e, -3 + 2 * .e), 50L, 4L, 1e-3, 1e-4)
    expect_false(is.na(.r$iter))
    expect_gte(.r$iter, 53L)
    expect_true(all(.r$drift <= .r$tol))
    expect_true(all(.r$tol > 1e-3 * c(5, 3)))
  })

  test_that("one trending coordinate holds the stop for all", {
    .x <- cbind(1, seq(0, 1, length.out = 300))
    expect_true(is.na(saemAutoStopTest_(.x, 50L, 1L, 1e-3, 1e-4)$iter))
  })

  test_that("saemControl() validates the auto-stop options", {
    expect_error(saemControl(autoStop = NA))
    expect_error(saemControl(stopWindow = 3))
    expect_error(saemControl(stopConsecutive = 0))
    expect_error(saemControl(stopRtol = -1))
    .c <- saemControl(autoStop = TRUE, stopWindow = 30)
    expect_true(.c$autoStop)
    expect_identical(.c$stopWindow, 30L)
    expect_error(saemControl(saCovAutoStop = NA))
    expect_true(saemControl()$saCovAutoStop)
    expect_false(saemControl(saCovAutoStop = FALSE)$saCovAutoStop)
  })

  .asMod <- function() {
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

  test_that("a simple PK fit stops before nEm and lands near the full run", {
    skip_on_cran()
    .ctl <- function(...) {
      saemControl(print = 0, nBurn = 100, nEm = 300, seed = 42L, calcTables = FALSE, ...)
    }
    .full <- suppressMessages(nlmixr2(.asMod(), theo_sd, est = "saem", control = .ctl(covMethod = "")))
    .stop <- suppressMessages(nlmixr2(.asMod(), theo_sd, est = "saem", control = .ctl(covMethod = "", autoStop = TRUE)))
    .as <- .stop$saemAutoStop
    expect_false(is.na(.as$iter))
    expect_lt(.as$iter, 400L)
    # the window only sees EM iterations
    expect_gte(.as$iter, 100L + 50L + 4L - 1L)
    expect_equal(.as$niter, 400L)
    expect_equal(names(.as$drift), c("tka", "tcl", "tv", "V(eta.ka)", "V(eta.cl)", "V(eta.v)", "add.sd"))
    expect_true(all(.as$drift <= .as$tol))
    # every per-iteration history ends at the stop
    expect_equal(nrow(.stop$parHistData), .as$iter)
    expect_equal(nrow(.stop$mcmcAccept), .as$iter)
    expect_equal(dim(.stop$phiM)[3], .as$iter)
    expect_equal(fixef(.stop), fixef(.full), tolerance = 0.02)
    expect_null(.full$saemAutoStop)
  })

  test_that("the SA covariance phase still runs after an early stop", {
    skip_on_cran()
    .f <- suppressMessages(nlmixr2(
      .asMod(),
      theo_sd,
      est = "saem",
      control = saemControl(
        print = 0,
        nBurn = 100,
        nEm = 300,
        seed = 42L,
        calcTables = FALSE,
        covMethod = "sa",
        nSaCov = 100L,
        autoStop = TRUE
      )
    ))
    expect_lt(.f$saemAutoStop$iter, 400L)
    expect_equal(.f$covMethod, "sa")
    expect_true(all(is.finite(diag(.f$cov))))
    expect_equal(nrow(.f$parHistData), .f$saemAutoStop$iter)
    .sc <- .f$saemSaCovStop
    expect_equal(.sc$nSaCov, 100L)
    expect_true(is.na(.sc$iter) || .sc$iter <= 100L)
  })

  test_that("the SA covariance phase stops on the same window and keeps the SEs", {
    skip_on_cran()
    .ctl <- function(...) {
      saemControl(
        print = 0,
        nBurn = 100,
        nEm = 300,
        seed = 42L,
        calcTables = FALSE,
        covMethod = "sa",
        nSaCov = 500L,
        ...
      )
    }
    .on <- suppressMessages(nlmixr2(.asMod(), theo_sd, est = "saem", control = .ctl()))
    .off <- suppressMessages(nlmixr2(.asMod(), theo_sd, est = "saem", control = .ctl(saCovAutoStop = FALSE)))
    .sc <- .on$saemSaCovStop
    expect_false(is.na(.sc$iter))
    expect_lt(.sc$iter, 500L)
    # a 10% burn-in, then a full window before the first stop
    expect_gte(.sc$iter, 50L + 50L + 4L - 1L)
    expect_true(all(.sc$drift <= .sc$tol))
    expect_equal(.on$covMethod, "sa")
    expect_equal(sqrt(diag(.on$cov)), sqrt(diag(.off$cov)), tolerance = 0.05)
    expect_equal(fixef(.on), fixef(.off))
    expect_null(.off$saemSaCovStop)
  })
})
