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

  test_that("a stationary point is left alone at a cost of 2n + 4 evaluations", {
    .n <- 0L
    .fn <- function(x) {
      .n <<- .n + 1L
      .valley(x)
    }
    .at <- list(par = c(1, 1), fval = 0, feval = 40L, ierr = 0L)
    .r <- .bobyqaStationary(.fn, .lo, .hi, .ctl, .at, tol = 0.01)
    expect_identical(.r$nStationaryRestart, 0L)
    expect_identical(.r$par, c(1, 1))
    expect_identical(.n, 8L)
    expect_identical(.r$feval, 48L)
  })

  test_that("a decrease below tol leaves the exit point alone", {
    .flat <- function(x) 1e-3 * sum((x - 1)^2)
    .at <- list(par = c(0, 0), fval = .flat(c(0, 0)), feval = 40L, ierr = 0L)
    .r <- .bobyqaStationary(.flat, .lo, .hi, .ctl, .at, tol = 0.01)
    expect_identical(.r$nStationaryRestart, 0L)
    expect_identical(.r$par, c(0, 0))
    expect_identical(.r$fval, .at$fval)
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

  test_that("the probe is held and kept out of the history", {
    .log <- character(0)
    .rec <- TRUE
    .record <- function(x) {
      .old <- .rec
      .rec <<- x
      .log <<- c(.log, paste0("record", x))
      .old
    }
    .hold <- function(x) {
      .log <<- c(.log, paste0("hold", x))
      FALSE
    }
    .fn <- function(x) {
      # every probe runs held and unrecorded; the restart runs released and recorded
      .log <<- c(.log, if (.rec) "fnRec" else "fnProbe")
      .valley(x)
    }
    .stop <- list(par = c(-0.5, 0.25), fval = .valley(c(-0.5, 0.25)), feval = 40L, ierr = 0L)
    .r <- .bobyqaStationary(.fn, .lo, .hi, .ctl, .stop, tol = 0.01, record = .record, hold = .hold)
    expect_gte(.r$nStationaryRestart, 1L)
    expect_true(.rec)
    .holds <- .log[startsWith(.log, "hold")]
    expect_identical(.holds[1], "holdTRUE")
    expect_identical(.holds[length(.holds)], "holdFALSE")
    # a probe's evaluations sit between a hold and its release
    .on <- 0L
    for (.e in .log) {
      if (.e == "holdTRUE") {
        .on <- 1L
      }
      if (.e == "holdFALSE") {
        .on <- 0L
      }
      if (.e == "fnProbe") {
        expect_identical(.on, 1L)
      }
      if (.e == "fnRec") expect_identical(.on, 0L)
    }
    expect_true(any(.log == "fnRec"))
    # a caller that had recording off keeps it off
    .rec <- FALSE
    .log <- character(0)
    .bobyqaStationary(
      .fn,
      .lo,
      .hi,
      .ctl,
      list(par = c(1, 1), fval = 0, feval = 40L, ierr = 0L),
      tol = 0.01,
      record = .record,
      hold = .hold
    )
    expect_false(.rec)
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

  test_that("a non-finite probe on a bound stops the check", {
    .nan <- function(x) if (x[1] > -5) NaN else -x[1] + x[2]^2
    .at <- list(par = c(-5, 0), fval = 5, feval = 10L, ierr = 0L)
    .r <- .bobyqaStationary(.nan, .lo, .hi, .ctl, .at, tol = 0.01)
    expect_identical(.r$nStationaryRestart, 0L)
    expect_identical(.r$par, c(-5, 0))
  })

  test_that("a restart that uses its whole budget stays within maxfun", {
    local_mocked_bindings(
      bobyqa = function(par, fn, control, ...) {
        list(par = c(1, 1), fval = 0, feval = control$maxfun, ierr = 1L)
      },
      .package = "minqa"
    )
    .stop <- list(par = c(-0.5, 0.25), fval = .valley(c(-0.5, 0.25)), feval = 40L, ierr = 0L)
    .r <- .bobyqaStationary(.valley, .lo, .hi, replace(.ctl, "maxfun", 200), .stop, tol = 0.01)
    expect_identical(.r$nStationaryRestart, 1L)
    expect_equal(.r$feval, 200)
  })

  test_that("restarts stop at maxRestart", {
    local_mocked_bindings(
      bobyqa = function(par, fn, ...) {
        .p <- par - c(0.05, 0)
        list(par = .p, fval = fn(.p), feval = 5L, ierr = 0L)
      },
      .package = "minqa"
    )
    .lin <- function(x) x[1]
    .at <- list(par = c(0, 0), fval = 0, feval = 10L, ierr = 0L)
    .r <- .bobyqaStationary(.lin, .lo, .hi, .ctl, .at, tol = 0.01)
    expect_identical(.r$nStationaryRestart, 3L)
  })

  test_that("infinite bounds and no budget cap still restart", {
    .stop <- list(par = c(-0.5, 0.25), fval = .valley(c(-0.5, 0.25)), feval = 40L, ierr = 0L)
    .r <- .bobyqaStationary(
      .valley,
      c(-Inf, -Inf),
      c(Inf, Inf),
      .ctl[names(.ctl) != "maxfun"],
      .stop,
      tol = 0.01
    )
    expect_gte(.r$nStationaryRestart, 1L)
    expect_lt(.r$fval, 1e-6)
  })

  test_that(".bobyqa() runs the check only when bobyqaStationary is set", {
    .base <- list(rhobeg = 0.2, rhoend = 1e-6, maxfun = 5000, sigdig = 3)
    .on <- .bobyqa(c(-1.2, 1), .valley, lower = .lo, upper = .hi, control = c(.base, bobyqaStationary = TRUE))
    expect_false(is.null(.on$nStationaryRestart))
    expect_lt(.on$fval, 1e-6)
    .off <- .bobyqa(c(-1.2, 1), .valley, lower = .lo, upper = .hi, control = c(.base, bobyqaStationary = FALSE))
    expect_null(.off$nStationaryRestart)
  })

  test_that("a fit the probe finds converged is the fit without the probe", {
    .one <- function() {
      ini({
        tka <- 0.45
        tcl <- 1.0
        tv <- 3.45
        add.sd <- 0.7
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .fit <- function(st) {
      suppressWarnings(nlmixr2(
        .one,
        nlmixr2data::theo_sd,
        est = "focei",
        control = foceiControl(print = 0, calcTables = FALSE, bobyqaStationary = st)
      ))
    }
    .on <- .fit(TRUE)
    .off <- .fit(FALSE)
    .o <- .on$env$optReturn
    expect_identical(.o$nStationaryRestart, 0L)
    # the probe ran (2n + 4 evaluations) and left the inner state as the search did
    expect_equal(.o$feval - .off$env$optReturn$feval, 2 * length(.o$par) + 4)
    expect_identical(.on$objf, .off$objf)
    expect_identical(.on$theta, .off$theta)
    expect_identical(.on$eta, .off$eta)
    expect_identical(.on$cov, .off$cov)
  })
})
