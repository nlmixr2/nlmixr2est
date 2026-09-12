test_that("nlminb consumes fast curvature and falls back only for Hessian errors", {
  f <- function(x) 100*(x[2]-x[1]^2)^2+(1-x[1])^2
  g <- function(x) c(-400*x[1]*(x[2]-x[1]^2)-2*(1-x[1]), 200*(x[2]-x[1]^2))
  for (mode in c("fast", "slow", "missing", "initialFail", "lateFail")) {
    calls <- 0L; messages <- character()
    hessian <- function(x) {
      calls <<- calls+1L
      if (mode == "initialFail" || (mode == "lateFail" && calls > 1L)) stop("no curvature")
      matrix(c(1200*x[1]^2-400*x[2]+2, -400*x[1], -400*x[1], 200), 2)
    }
    result <- withCallingHandlers(.nlminb(c(-1.2, 1), f, g,
      control = list(fast = mode != "slow", maxOuterIterations = 1000L,
                     hessian = if (mode == "missing") NULL else hessian)),
      warning = function(w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") })
    expect_equal(result$x, c(1, 1), tolerance = 1e-5)
    expect_equal(result$hessianEvaluations, calls)
    expect_identical(result$hessianFallback, mode %in% c("initialFail", "lateFail"))
    if (mode %in% c("slow", "missing")) expect_identical(calls, 0L)
    else expect_gt(calls, 0L)
    if (mode %in% c("initialFail", "lateFail")) {
      expect_length(messages, 1L)
      expect_match(messages, "restarting gradient-only nlminb")
    } else expect_length(messages, 0L)
  }
  expect_error(.nlminb(c(0, 0), function(x) stop("objective failure"), g,
    control = list(fast = TRUE, hessian = function(x) diag(2))), "objective failure")
})
