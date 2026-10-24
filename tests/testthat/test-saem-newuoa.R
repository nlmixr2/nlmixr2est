# .newuoa() is both the SAEM phi0 refinement and the outer optimizer for
# foceiControl(outerOpt = "newuoa"), so it must stop at a foceiControl() list's
# maxOuterIterations the way .uobyqa() does, and still honour a plain maxfun.
# Rosenbrock from (-1.2, 1) needs well over 50 evaluations unbounded at these
# settings (the exact count differs by platform's floating point: 160 on x86-64
# Linux, 141 on arm64 macOS), so a limit of 50 always binds; 50 is above minqa's
# 10 * length(par)^2 advice, so no warning is raised.

.rosenbrock <- function(x, ...) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2

test_that(".newuoa() stops at maxOuterIterations evaluations", {
  .ret <- .newuoa(c(-1.2, 1), .rosenbrock, control = list(maxOuterIterations = 50L))
  expect_identical(.ret$feval, 50L)
  expect_equal(.ret$convergence, 1)
})

test_that(".newuoa() and .uobyqa() apply maxOuterIterations alike", {
  .ctl <- list(maxOuterIterations = 50L)
  expect_identical(
    .newuoa(c(-1.2, 1), .rosenbrock, control = .ctl)$feval,
    .uobyqa(c(-1.2, 1), .rosenbrock, control = .ctl)$feval
  )
})

test_that(".newuoa() honours maxfun when no maxOuterIterations is given", {
  .ret <- .newuoa(c(-1.2, 1), .rosenbrock, control = list(maxfun = 50L))
  expect_identical(.ret$feval, 50L)
  expect_equal(.ret$convergence, 1)
})

test_that(".newuoa() without a limit runs past 50 evaluations to the minimum", {
  .ret <- .newuoa(c(-1.2, 1), .rosenbrock)
  expect_gt(.ret$feval, 50L)
  expect_equal(.ret$convergence, 0)
  expect_equal(.ret$x, c(1, 1), tolerance = 1e-3)
})
