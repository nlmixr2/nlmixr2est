test_that(".covPrecursorCheck() takes NULL or known sources, in order", {
  expect_identical(.covPrecursorCheck(NULL), character(0))
  expect_identical(.covPrecursorCheck(character(0)), character(0))
  expect_identical(.covPrecursorCheck("analytic"), "analytic")
  expect_identical(.covPrecursorCheck(c("analytic", "fd")), c("analytic", "fd"))
  expect_error(.covPrecursorCheck("saem"), "unknown source\\(s\\) \"saem\"; use some of \"fd\", \"analytic\", or NULL for none")
  expect_error(.covPrecursorCheck(c("fd", "fd")), "names a source more than once")
  expect_error(.covPrecursorCheck(NA_character_), "must be NULL or some of")
  expect_error(.covPrecursorCheck(1), "must be NULL or some of")
})

test_that("foceiControl() and rsControl() carry covPrecursor", {
  expect_identical(foceiControl()$covPrecursor, c("fd", "analytic"))
  expect_identical(foceiControl(covPrecursor = NULL)$covPrecursor, character(0))
  expect_identical(foceiControl(covPrecursor = "analytic")$covPrecursor, "analytic")
  expect_error(foceiControl(covPrecursor = "nlme"), "unknown source")
  # rsControl(): left out keeps the fit's value, NULL is none
  expect_false("covPrecursor" %in% names(rsControl()))
  expect_identical(rsControl(covPrecursor = NULL)$covPrecursor, character(0))
  expect_identical(rsControl(covPrecursor = "fd")$covPrecursor, "fd")
})

test_that("covPrecursor deparses and round-trips", {
  for (.v in list(character(0), "analytic", c("analytic", "fd"))) {
    .back <- eval(rxode2::rxUiDeparse(foceiControl(covPrecursor = .v), "ctl")[[3]])
    expect_identical(.back$covPrecursor, .v)
  }
  # the default is not written out
  expect_false(grepl("covPrecursor", deparse1(rxode2::rxUiDeparse(foceiControl(), "ctl")), fixed = TRUE))
  .back <- eval(rxode2::rxUiDeparse(rsControl(covPrecursor = NULL), "ctl")[[3]])
  expect_identical(.back$covPrecursor, character(0))
})

test_that("covPrecursor is in the covariance store's key", {
  expect_true("covPrecursor" %in% .covStoreKeyFields)
  expect_true(.covStoreRefitOk(list(covPrecursor = "fd")))
})

test_that(".covPrecursorLine() describes how a precursor served", {
  expect_null(.covPrecursorLine(NULL))
  expect_identical(.covPrecursorLine(list(source = "fd")), "from the \"fd\" precursor (seeded steps)")
})

test_that(".covPrecursorFd() reads a full R stored for other settings at the same estimates", {
  .nm <- c("a", "b")
  .R <- matrix(c(2, 0.1, 0.1, 3), 2, dimnames = list(.nm, .nm))
  .key <- list(handoff = list(theta = 1), settings = list(hessEps = 1), versions = c(x = "1"))
  .env <- new.env(parent = emptyenv())
  # the entry for the very same settings is the exact store, not a hint
  .env$covStore <- list(list(key = .key, full = list(R = .R)))
  expect_null(.covPrecursorFd(.env, .key, .nm))
  .other <- .key
  .other$settings$hessEps <- 2
  .env$covStore <- list(list(key = .other, full = list(R = .R)))
  expect_identical(.covPrecursorFd(.env, .key, .nm), .R)
  # another hand-off, other versions or other parameters are not this point
  .moved <- .other
  .moved$handoff$theta <- 2
  .env$covStore <- list(list(key = .moved, full = list(R = .R)))
  expect_null(.covPrecursorFd(.env, .key, .nm))
  .old <- .other
  .old$versions <- c(x = "0")
  .env$covStore <- list(list(key = .old, full = list(R = .R)))
  expect_null(.covPrecursorFd(.env, .key, .nm))
  .env$covStore <- list(list(key = .other, full = list(R = .R)))
  expect_null(.covPrecursorFd(.env, .key, c("a", "c")))
  expect_null(.covPrecursorFd(.env, NULL, .nm))
})

test_that(".covPrecursorAnalytic() inverts an analytic covariance over the full parameters", {
  .nm <- c("a", "b")
  .cov <- matrix(c(0.5, 0.1, 0.1, 0.4), 2, dimnames = list(.nm, .nm))
  .env <- new.env(parent = emptyenv())
  .env$cov <- .cov
  .env$covMethod <- "analytic (full)"
  expect_equal(.covPrecursorAnalytic(.env, .nm), solve(.cov))
  # from covList as well, and only an analytic one over these parameters
  .env$covMethod <- "r,s (full)"
  expect_null(.covPrecursorAnalytic(.env, .nm))
  .env$covList <- list(analytic = .cov)
  expect_equal(.covPrecursorAnalytic(.env, .nm), solve(.cov))
  expect_null(.covPrecursorAnalytic(.env, c("a", "c")))
  .env$covList <- list(analytic = matrix(0, 2, 2, dimnames = list(.nm, .nm)))
  expect_null(.covPrecursorAnalytic(.env, .nm))
})

nmTest({
  .pcModel <- function() {
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
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  test_that("setCov() starts the full stage from the analytic covariance the fit holds", {
    skip_on_cran()
    .ctl <- foceiControl(print = 0, calcTables = FALSE, covMethod = "analytic")
    .none <- suppressMessages(.nlmixr(.pcModel, theo_sd, "focei", .ctl))
    .pre <- suppressMessages(.nlmixr(.pcModel, theo_sd, "focei", .ctl))
    suppressMessages(suppressWarnings(setCov(.none, "r,s (full)", rsControl(covPrecursor = NULL))))
    suppressMessages(suppressWarnings(setCov(.pre, "r,s (full)", rsControl(covPrecursor = "analytic"))))
    expect_identical(.none$covMethod, "r,s (full)")
    expect_identical(.pre$covMethod, "r,s (full)")
    expect_null(.none$env$covPrecursorUsed[["r,s (full)"]])
    expect_identical(.pre$env$covPrecursorUsed[["r,s (full)"]]$source, "analytic")
    # seeded steps pass the same acceptance test: the same covariance to within the
    # finite-difference error
    expect_equal(sqrt(diag(.pre$cov)), sqrt(diag(.none$cov)), tolerance = 0.01)
    expect_output(print(.pre), "from the \"analytic\" precursor (seeded steps", fixed = TRUE)
  })
})
