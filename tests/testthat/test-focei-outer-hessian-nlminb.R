test_that("nlminb consumes fast curvature and falls back only for Hessian errors", {
  f <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
  g <- function(x) c(-400 * x[1] * (x[2] - x[1]^2) - 2 * (1 - x[1]), 200 * (x[2] - x[1]^2))
  for (mode in c("fast", "slow", "missing", "initialFail", "lateFail")) {
    calls <- 0L
    messages <- character()
    hessian <- function(x) {
      calls <<- calls + 1L
      if (mode == "initialFail" || (mode == "lateFail" && calls > 1L)) {
        stop("no curvature")
      }
      matrix(c(1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1], -400 * x[1], 200), 2)
    }
    result <- withCallingHandlers(
      .nlminb(
        c(-1.2, 1),
        f,
        g,
        control = list(
          fast = mode != "slow",
          maxOuterIterations = 1000L,
          hessian = if (mode == "missing") NULL else hessian
        )
      ),
      warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_equal(result$x, c(1, 1), tolerance = 1e-5)
    expect_equal(result$hessianEvaluations, calls)
    expect_identical(result$hessianFallback, mode %in% c("initialFail", "lateFail"))
    if (mode %in% c("slow", "missing")) {
      expect_identical(calls, 0L)
    } else {
      expect_gt(calls, 0L)
    }
    if (mode %in% c("initialFail", "lateFail")) {
      expect_length(messages, 1L)
      expect_match(messages, "restarting gradient-only nlminb")
    } else {
      expect_length(messages, 0L)
    }
  }
  expect_error(
    .nlminb(
      c(0, 0),
      function(x) stop("objective failure"),
      g,
      control = list(fast = TRUE, hessian = function(x) diag(2))
    ),
    "objective failure"
  )
})

test_that('nlminb curvature is opt-in for normalized fast controls', {
  objective <- function(x) sum(x^2)
  gradient <- function(x) 2 * x
  for (mode in c('auto', 'bfgs', 'analytic', 'fd')) {
    calls <- 0L
    hessian <- function(x) {
      calls <<- calls + 1L
      diag(2, length(x))
    }
    result <- .nlminb(
      c(2, 3),
      objective,
      gradient,
      control = list(fast = TRUE, outerTrustHessian = mode, hessian = hessian, maxOuterIterations = 100L)
    )
    expect_equal(result$x, c(0, 0), tolerance = 1e-7)
    if (mode == 'analytic') {
      expect_gt(calls, 0L)
    } else {
      expect_identical(calls, 0L)
    }
    if (mode == 'fd') expect_gt(result$hessianEvaluations, 0L)
  }
})

test_that('nlminb FD curvature does not need the analytical Hessian', {
  for (fast in c(FALSE, TRUE)) {
    result <- .nlminb(
      c(2, 3),
      function(x) sum(x^2),
      function(x) 2 * x,
      control = list(
        fast = fast,
        outerTrustHessian = 'fd',
        hessian = function(x) stop('analytical Hessian must not be called'),
        maxOuterIterations = 100L
      )
    )
    expect_equal(result$x, c(0, 0), tolerance = 1e-7)
    expect_gt(result$hessianEvaluations, 0L)
    expect_false(result$hessianFallback)
  }
  expect_warning(
    result <- .nlminb(
      0,
      function(x) (x - 2e-7)^2,
      function(x) 2 * (x - 2e-7),
      lower = -1e-6,
      upper = 1e-6,
      control = list(outerTrustHessian = 'fd', maxOuterIterations = 100L)
    ),
    'restarting gradient-only nlminb'
  )
  expect_true(result$hessianFallback)
  expect_true(all(is.finite(result$x)))
})

test_that('nlminb FD captures the base gradient before moving the objective', {
  point <- NULL
  fn <- function(x) {
    point <<- x
    sum(x^2)
  }
  gr <- function(x) {
    if (!identical(point, x)) {
      stop('gradient point differs from objective point')
    }
    2 * x
  }
  result <- .nlminb(c(2, 3), fn, gr, control = list(outerTrustHessian = 'fd', maxOuterIterations = 100L))
  expect_equal(result$x, c(0, 0), tolerance = 1e-7)
  expect_gt(result$hessianEvaluations, 0L)
  expect_false(result$hessianFallback)
})

test_that('nlminb recovers when an FD probe has a nonfinite gradient', {
  gr <- function(x) if (x[1] > 2) c(NA_real_, NA_real_) else 2 * x
  expect_warning(
    result <- .nlminb(
      c(2, 3),
      function(x) sum(x^2),
      gr,
      control = list(outerTrustHessian = 'fd', maxOuterIterations = 100L)
    ),
    'restarting gradient-only nlminb'
  )
  expect_equal(result$x, c(0, 0), tolerance = 1e-7)
  expect_true(result$hessianFallback)
})
