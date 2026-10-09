# The finite-difference covariance options foceiControl() and rsControl() share.

test_that("foceiControl() and rsControl() take the same finite-difference covariance options", {
  # gillStepCov is the factor the Gill search grows its step by
  expect_error(foceiControl(gillStepCov = 0.5), "gillStepCov")
  expect_error(rsControl(gillStepCov = 0.5), "gillStepCov")
  expect_identical(foceiControl(gillStepCov = 1)$gillStepCov, 1)
  expect_identical(rsControl(gillStepCov = 1)$gillStepCov, 1)
  # covSmall is a single finite threshold
  expect_error(foceiControl(covSmall = Inf), "covSmall")
  expect_error(rsControl(covSmall = Inf), "covSmall")
  expect_error(foceiControl(covSmall = c(1e-5, 1e-6)), "covSmall")
  expect_error(rsControl(covSmall = c(1e-5, 1e-6)), "covSmall")
  # the flags are logical or 0/1, which foceiControl() stores
  expect_identical(
    unclass(rsControl(covGillF = 0L, rmatNorm = 1, smatNorm = TRUE)),
    list(covGillF = 0L, rmatNorm = 1, smatNorm = TRUE)
  )
  expect_identical(
    foceiControl(covGillF = 0L, rmatNorm = 1, smatNorm = TRUE)[c("covGillF", "rmatNorm", "smatNorm")],
    list(covGillF = 0L, rmatNorm = 1L, smatNorm = 1L)
  )
  expect_error(rsControl(rmatNorm = 2), "rmatNorm")
  expect_error(foceiControl(rmatNorm = 2), "rmatNorm")
})

test_that("every Gill step factor is a finite number of at least 1", {
  # gillStep (the fit's gradient search), gillStepCov and gillStepCovLlik (the
  # covariance step) all grow the Gill (1983) step by multiplying by the factor
  # and shrink it by dividing
  for (.n in c("gillStep", "gillStepCov", "gillStepCovLlik")) {
    .below <- paste0("Assertion on '", .n, "' failed: Element 1 is not >= 1.")
    .inf <- paste0("Assertion on '", .n, "' failed: Must be finite.")
    expect_error(do.call(foceiControl, setNames(list(0.5), .n)), .below, fixed = TRUE)
    expect_error(do.call(foceiControl, setNames(list(0), .n)), .below, fixed = TRUE)
    expect_error(do.call(foceiControl, setNames(list(Inf), .n)), .inf, fixed = TRUE)
    expect_identical(do.call(foceiControl, setNames(list(1), .n))[[.n]], 1)
  }
  expect_error(rsControl(gillStepCov = Inf), "Assertion on 'gillStepCov' failed: Must be finite.", fixed = TRUE)
  # nlmixr2Gill83() makes the checks foceiControl() makes
  .f <- function(x) sum(x^2)
  expect_error(
    nlmixr2Gill83(.f, c(1, 2), gillStep = 0.5),
    "Assertion on 'gillStep' failed: Element 1 is not >= 1.",
    fixed = TRUE
  )
  expect_error(
    nlmixr2Gill83(.f, c(1, 2), gillStep = Inf),
    "Assertion on 'gillStep' failed: Must be finite.",
    fixed = TRUE
  )
})

test_that("rsControl(covFallback=) is the request's own list, never the fit's", {
  expect_null(rsControl()$covFallback)
  expect_identical(rsControl(covFallback = list(r = "s"))$covFallback, list(r = "s"))
  expect_error(rsControl(covFallback = list(r = "vi")), "cannot fall back to \"vi\"")
  .fit <- new.env(parent = emptyenv())
  .fit$foceiControl <- foceiControl()
  expect_null(setCovOptions(rsControl(), .fit)$covFallback)
  expect_identical(setCovOptions(rsControl(covFallback = list(r = "s")), .fit)$covFallback, list(r = "s"))
  .d <- rxode2::rxUiDeparse(rsControl(covFallback = list(r = "s")), "ctl")
  expect_identical(eval(.d[[3]])$covFallback, list(r = "s"))
})

nmTest({
  test_that("setCov() falls back only as rsControl(covFallback=) lists", {
    skip_on_cran()
    .m <- function() {
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
    .fit <- .nlmixr(.m, theo_sd, "focei", foceiControl(print = 0, calcTables = FALSE, covMethod = "r,s", covFull = FALSE))
    # an "r" request whose refit can only give "s"
    .s <- suppressMessages(.setCovRefit(.fit, covMethod = "s", covFull = FALSE))
    expect_identical(.s$covMethod, "s")
    .args <- new.env(parent = emptyenv())
    local_mocked_bindings(.setCovRefit = function(obj, ...) {
      .args$fallback <- list(...)$covFallback
      .s
    })
    expect_error(suppressMessages(setCov(.fit, "r")), "\"r\" could not be computed")
    # the refit was given no fallback
    expect_identical(.args$fallback, list())
    expect_warning(
      suppressMessages(setCov(.fit, "r", control = rsControl(covFallback = list(r = "s")))),
      "\"s\" covariance installed instead of the requested \"r\""
    )
    expect_identical(.args$fallback, list(r = "s"))
    expect_identical(.fit$covMethod, "s")
    expect_identical(.fit$cov, .s$cov)
  })
})
