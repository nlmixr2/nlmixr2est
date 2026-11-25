test_that("the stop-file poll removes the file and stays requested until reset", {
  .sf <- tempfile("stop")
  on.exit({
    .nmStopFileSet("")
    unlink(.sf)
  })
  .nmStopFileSet(.sf)
  expect_false(.nmStopFilePoll())
  file.create(.sf)
  expect_true(.nmStopFilePoll())
  expect_false(file.exists(.sf))
  expect_true(.nmStopFileRequested())
  expect_true(.nmStopFileStopped())
  .nmStopFileSet("")
  expect_false(.nmStopFileRequested())
  expect_false(.nmStopFileStopped())
  file.create(.sf)
  expect_false(.nmStopFilePoll())
})

test_that("a nested estimation keeps the outer stop file", {
  .sf <- tempfile("stop")
  on.exit(.nmStopFileSet(""))
  expect_false(.nmStopFileStart(NULL))
  expect_false(.nmStopFileStart(list(stopFile = NULL)))
  expect_true(.nmStopFileStart(list(stopFile = .sf)))
  expect_false(.nmStopFileStart(list(stopFile = tempfile("other"))))
  expect_equal(.nmStopFileGet(), normalizePath(.sf, mustWork = FALSE))
  expect_false(.nmStopFileEnd())
  expect_equal(.nmStopFileGet(), "")
})

test_that("foceiControl(stopFile=) is validated", {
  expect_null(foceiControl()$stopFile)
  expect_equal(foceiControl(stopFile = "stop")$stopFile, "stop")
  expect_error(foceiControl(stopFile = ""))
  expect_error(foceiControl(stopFile = 1))
})

skip_on_cran()

.stopMod <- function() {
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

for (.oo in c("bobyqa", "nlminb", "L-BFGS-B", "lbfgsb3c")) {
  test_that(paste0("a stop file ends a FOCEi fit early (outerOpt=", .oo, ")"), {
    .sf <- tempfile("stop")
    on.exit(unlink(.sf))
    file.create(.sf)
    .fit <- suppressMessages(nlmixr(
      .stopMod,
      nlmixr2data::theo_sd,
      "focei",
      control = foceiControl(print = 0, stopFile = .sf, outerOpt = .oo, covMethod = "")
    ))
    expect_false(file.exists(.sf))
    expect_equal(.fit$message, "stopped early by user request")
    expect_true("stopped early by user request" %in% .fit$runInfo)
    # stopped at the first evaluation: the initial estimates
    expect_equal(.fit$theta[["tka"]], 0.45, tolerance = 1e-6)
    expect_equal(.nmStopFileGet(), "")
    expect_false(.nmStopFileRequested())
  })
}

test_that("a stop part way through keeps the lowest evaluation and a re-run is not stopped", {
  .sf <- tempfile("stop")
  on.exit(unlink(.sf))
  .env <- new.env(parent = emptyenv())
  .env$n <- 0L
  .opt <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
    .fn <- function(x) {
      .ret <- fn(x)
      .env$n <- .env$n + 1L
      if (.env$n == 5L) {
        file.create(.sf)
      }
      .ret
    }
    .bobyqa(par, .fn, gr, lower = lower, upper = upper, control = control, ...)
  }
  .fit <- suppressMessages(nlmixr(
    .stopMod,
    nlmixr2data::theo_sd,
    "focei",
    control = foceiControl(print = 0, stopFile = .sf, outerOpt = .opt, covMethod = "")
  ))
  expect_false(file.exists(.sf))
  expect_equal(.env$n, 5L)
  expect_equal(.fit$message, "stopped early by user request")
  expect_true("stopped early by user request" %in% .fit$runInfo)
  .ph <- .fit$parHistData
  .obj <- .ph$objf[.ph$type == "Unscaled"]
  expect_length(.obj, 6L)
  expect_equal(.fit$objf, min(.obj), tolerance = 1e-3)

  .fit2 <- suppressMessages(nlmixr(
    .stopMod,
    nlmixr2data::theo_sd,
    "focei",
    control = foceiControl(print = 0, stopFile = .sf, covMethod = "")
  ))
  expect_false(.fit2$message == "stopped early by user request")
  expect_false("stopped early by user request" %in% .fit2$runInfo)
})

test_that("a stop file ends a SAEM fit at the current iteration", {
  expect_null(saemControl()$stopFile)
  expect_error(saemControl(stopFile = ""))
  .sf <- tempfile("stop")
  on.exit(unlink(.sf))
  file.create(.sf)
  .fit <- suppressMessages(nlmixr(
    .stopMod,
    nlmixr2data::theo_sd,
    "saem",
    control = saemControl(print = 0, stopFile = .sf, nBurn = 20, nEm = 30)
  ))
  expect_false(file.exists(.sf))
  expect_true("stopped early by user request" %in% .fit$runInfo)
  expect_equal(nrow(.fit$parHist), 1L)
  expect_true(is.finite(.fit$objf))
  expect_equal(.nmStopFileGet(), "")

  .fit2 <- suppressMessages(nlmixr(
    .stopMod,
    nlmixr2data::theo_sd,
    "saem",
    control = saemControl(print = 0, stopFile = .sf, nBurn = 20, nEm = 30)
  ))
  expect_false("stopped early by user request" %in% .fit2$runInfo)
  expect_equal(nrow(.fit2$parHist), 50L)
})
