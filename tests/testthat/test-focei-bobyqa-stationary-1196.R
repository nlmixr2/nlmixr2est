nmTest({
  # bobyqa shrank its trust region in a narrow valley and exited normally short of
  # the minimum (#1196); a gradient probe now restarts it.
  .valley <- function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
  .ctl <- list(npt = 5, rhobeg = 0.2, rhoend = 1e-6, iprint = 0L, maxfun = 5000)
  .lo <- c(-5, -5)
  .hi <- c(5, 5)

  test_that("bobyqaStationary is on by default and round-trips", {
    expect_true(foceiControl()$bobyqaStationary)
    expect_false(foceiControl(bobyqaStationary = FALSE)$bobyqaStationary)
    expect_error(foceiControl(bobyqaStationary = NA))
    .dep <- rxUiDeparse(foceiControl(bobyqaStationary = FALSE), "ctl")
    expect_false(eval(.dep[[3]])$bobyqaStationary)
  })

  test_that("a normal exit short of the minimum is restarted", {
    .stop <- list(par = c(-0.5, 0.25), fval = .valley(c(-0.5, 0.25)), feval = 40L, ierr = 0L)
    .r <- .bobyqaStationary(.valley, .lo, .hi, .ctl, .stop, tol = 0.01)
    expect_gte(.r$nStationaryRestart, 1L)
    expect_lt(.r$fval, 1e-6)
    expect_equal(.r$par, c(1, 1), tolerance = 1e-3)
    expect_gt(.r$feval, 40L)
  })

  test_that("a stationary point is left alone at a cost of 2n + 5 evaluations", {
    .n <- 0L
    .fn <- function(x) {
      .n <<- .n + 1L
      .valley(x)
    }
    .at <- list(par = c(1, 1), fval = 0, feval = 40L, ierr = 0L)
    .r <- .bobyqaStationary(.fn, .lo, .hi, .ctl, .at, tol = 0.01)
    expect_identical(.r$nStationaryRestart, 0L)
    expect_identical(.r$par, c(1, 1))
    expect_identical(.n, 9L)
    expect_identical(.r$feval, 49L)
  })

  test_that("no check after a non-normal exit or without budget", {
    .n <- 0L
    .fn <- function(x) {
      .n <<- .n + 1L
      .valley(x)
    }
    .stop <- list(par = c(-0.5, 0.25), fval = .valley(c(-0.5, 0.25)), feval = 40L, ierr = 1L)
    expect_identical(.bobyqaStationary(.fn, .lo, .hi, .ctl, .stop, tol = 0.01)$par, .stop$par)
    .stop$ierr <- 0L
    expect_identical(.bobyqaStationary(.fn, .lo, .hi, replace(.ctl, "maxfun", 50), .stop, tol = 0.01)$par, .stop$par)
    expect_identical(.n, 0L)
  })

  test_that("a point on a bound sees an inward gradient and stays in the box", {
    .lin <- function(x) -x[1]
    .at <- list(par = c(-5, 0), fval = 5, feval = 10L, ierr = 0L)
    .r <- .bobyqaStationary(.lin, .lo, .hi, .ctl, .at, tol = 0.01)
    expect_gte(.r$nStationaryRestart, 1L)
    expect_equal(.r$par[1], 5)
    expect_true(all(.r$par >= .lo & .r$par <= .hi))
    # descent that would leave the box is not a reason to restart
    .out <- list(par = c(5, 0), fval = -5, feval = 10L, ierr = 0L)
    .r <- .bobyqaStationary(.lin, .lo, .hi, .ctl, .out, tol = 0.01)
    expect_identical(.r$nStationaryRestart, 0L)
    expect_identical(.r$par, c(5, 0))
  })

  test_that("a non-finite exit value is left alone", {
    .at <- list(par = c(-0.5, 0.25), fval = NaN, feval = 40L, ierr = 0L)
    .r <- .bobyqaStationary(.valley, .lo, .hi, .ctl, .at, tol = 0.01)
    expect_identical(.r$par, .at$par)
    expect_identical(.r$nStationaryRestart, 0L)
  })

  test_that("the probe point is kept when the restart does not beat it", {
    local_mocked_bindings(
      bobyqa = function(par, fn, ...) list(par = par + 1, fval = Inf, feval = 3L, ierr = 0L),
      .package = "minqa"
    )
    .stop <- list(par = c(-0.5, 0.25), fval = .valley(c(-0.5, 0.25)), feval = 40L, ierr = 0L)
    .r <- .bobyqaStationary(.valley, .lo, .hi, .ctl, .stop, tol = 0.01)
    expect_identical(.r$nStationaryRestart, 1L)
    expect_lt(.r$fval, .stop$fval - 0.01)
    expect_equal(.r$fval, .valley(.r$par))
    expect_true(all(abs(.r$par - .stop$par) < 0.03))
  })

  test_that(".bobyqa() runs the check only when bobyqaStationary is set", {
    .base <- list(rhobeg = 0.2, rhoend = 1e-6, maxfun = 5000, sigdig = 3)
    .on <- .bobyqa(c(-1.2, 1), .valley, lower = .lo, upper = .hi, control = c(.base, bobyqaStationary = TRUE))
    expect_false(is.null(.on$nStationaryRestart))
    expect_lt(.on$fval, 1e-6)
    .off <- .bobyqa(c(-1.2, 1), .valley, lower = .lo, upper = .hi, control = c(.base, bobyqaStationary = FALSE))
    expect_null(.off$nStationaryRestart)
  })
})
