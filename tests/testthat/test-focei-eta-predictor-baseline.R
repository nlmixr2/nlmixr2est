test_that('ETA extrapolation and its bounded fallback use the gradient-point baseline', {
  model <- function() {
    ini({ theta <- .2; slope <- .1; error <- .2; eta ~ .3 })
    model({ d/dt(dummy) <- 0; dummy(0) <- 1
            prediction <- exp(theta + eta) * dummy + slope * TIME
            prediction ~ add(error) })
  }
  data <- data.frame(ID = rep(1:2, each = 3), TIME = rep(1:3, 2), DV = rep(exp(.2) + .1 * (1:3), 2), AMT = 0, EVID = 0)
  for (reset in c(0, .5)) {
    values <- NULL
    optimizer <- function(par, fn, gr, ...) {
      fn(par)
      gr(par)
      a <- b <- par
      a[1] <- a[1] + .8
      b[1] <- b[1] - .2
      first <- fn(a)
      fn(b)
      values <<- c(first, fn(a))
      list(x = par, convergence = 0L, message = 'line-search regression')
    }
    # One inner step exposes the starting point; converged solves can hide it.
    fit <- suppressWarnings(suppressMessages(nlmixr2(
      model,
      data,
      'focei',
      foceiControl(
        fast = TRUE,
        outerOpt = optimizer,
        innerOpt = 'trust',
        maxOuterIterations = 10,
        maxInnerIterations = 1,
        mceta = -2,
        resetEtaP = reset,
        etaNudge = 0,
        etaNudge2 = 0,
        print = 0,
        covMethod = '',
        calcTables = FALSE
      )
    )))
    expect_length(values, 2)
    expect_equal(values[1], values[2], tolerance = 1e-10)
    expect_gt(fit$env$nAnalyticGradDirect, 0)
    expect_equal(fit$env$nFDGradFast, 0)
  }
})

test_that('an analytic gradient reuses the inner solution of its objective point', {
  model <- function() {
    ini({ theta <- .2; slope <- .1; error <- .2; eta ~ .3 })
    model({ d/dt(dummy) <- 0; dummy(0) <- 1
            prediction <- exp(theta + eta) * dummy + slope * TIME
            prediction ~ add(error) })
  }
  data <- data.frame(ID = rep(1:2, each = 3), TIME = rep(1:3, 2), DV = c(1.5, 1.7, 1.6, .8, 1, 1.1), AMT = 0, EVID = 0)
  observed <- list()
  for (withGradient in c(FALSE, TRUE)) {
    values <- NULL
    optimizer <- function(par, fn, gr, ...) {
      first <- fn(par)
      if (withGradient) {
        gr(par)
      }
      values <<- c(first, fn(par))
      list(x = par, convergence = 0L, message = 'gradient point regression')
    }
    fit <- suppressWarnings(suppressMessages(nlmixr2(
      model,
      data,
      'focei',
      foceiControl(
        fast = TRUE,
        outerOpt = optimizer,
        innerOpt = 'trust',
        maxOuterIterations = 10,
        maxInnerIterations = 1,
        mceta = -2,
        resetEtaP = 0,
        etaNudge = 0,
        etaNudge2 = 0,
        print = 0,
        covMethod = '',
        calcTables = FALSE
      )
    )))
    observed[[as.character(withGradient)]] <- values
    if (withGradient) {
      expect_gt(fit$env$nAnalyticGradDirect, 0)
      expect_equal(fit$env$nFDGradFast, 0)
    }
  }
  expect_equal(observed[["TRUE"]], observed[["FALSE"]], tolerance = 1e-10)
})

test_that('an exact analytic zero does not replace the selected outer optimizer', {
  model <- function() {
    ini({ theta <- .2; slope <- .1; error <- .2; eta ~ .3 })
    model({ d/dt(dummy) <- 0; dummy(0) <- 1
            prediction <- exp(theta + eta) * dummy + slope * Z
            prediction ~ add(error) })
  }
  z <- rep(c(-1, 1), 4)
  data <- data.frame(ID = rep(1:2, each = 4), TIME = rep(1:4, 2), Z = z, DV = exp(.2) + .1 * z, AMT = 0, EVID = 0)
  gradient <- NULL
  optimizer <- function(par, fn, gr, ...) {
    fn(par)
    gradient <<- gr(par)
    list(x = par, convergence = 0L, message = 'exact analytic zero regression')
  }
  fit <- suppressWarnings(suppressMessages(nlmixr2(
    model,
    data,
    'focei',
    foceiControl(
      fast = TRUE,
      outerOpt = optimizer,
      innerOpt = 'trust',
      maxOuterIterations = 10,
      etaNudge = 0,
      etaNudge2 = 0,
      zeroGradFirstReset = FALSE,
      print = 0,
      covMethod = '',
      calcTables = FALSE
    )
  )))
  expect_length(gradient, 4)
  expect_identical(gradient[2], 0)
  expect_gt(fit$env$nAnalyticGradDirect, 0)
  expect_equal(fit$env$nFDGradFast, 0)
  expect_identical(fit$message, 'exact analytic zero regression')
})

test_that('analytic point caching does not cross fits with the same theta', {
  model <- function() {
    ini({ level <- .2; error <- .2; eta ~ .3 })
    model({ d/dt(dummy) <- 0; dummy(0) <- 1
            prediction <- (level + eta) * dummy; prediction ~ add(error) })
  }
  for (value in c(1.6, 2)) {
    data <- data.frame(ID = 1, TIME = 1:3, DV = value, AMT = 0, EVID = 0)
    eta <- 3 * .3 * (value - .2) / (.2^2 + 3 * .3)
    fit <- suppressWarnings(suppressMessages(nlmixr2(
      model,
      data,
      'focei',
      foceiControl(
        fast = TRUE,
        maxOuterIterations = 0,
        maxInnerIterations = 0,
        etaMat = matrix(eta, 1, 1),
        resetEtaP = 0,
        print = 0,
        covMethod = '',
        calcTables = FALSE
      )
    )))
    expect_equal(fit$env$.gradDirectFirst[1], -2 * 3 * (value - .2) / (.2^2 + 3 * .3), tolerance = 1e-7)
    expect_gt(fit$env$nAnalyticGradDirect, 0)
    expect_equal(fit$env$nFDGradFast, 0)
  }
})

test_that('a flat rejected fast trial does not replace the outer optimizer', {
  model <- function() {
    ini({ level <- .2; error <- .2; eta ~ .3 })
    model({ d/dt(dummy) <- 0; dummy(0) <- 1
            prediction <- exp(level + eta) * dummy; prediction ~ add(error) })
  }
  data <- data.frame(ID = 1, TIME = 1:3, DV = 1.6, AMT = 0, EVID = 0)
  trialReturned <- FALSE
  values <- NULL
  optimizer <- function(par, fn, gr, ...) {
    initial <- fn(par)
    gr(par)
    invalid <- par
    invalid[1] <- 1000
    penalty <- fn(invalid)
    gr(invalid)
    trialReturned <<- TRUE
    values <<- c(initial, penalty, fn(par))
    list(x = par, convergence = 0L, message = 'rejected trial regression')
  }
  fit <- suppressWarnings(suppressMessages(nlmixr2(
    model,
    data,
    'focei',
    foceiControl(
      fast = TRUE,
      outerOpt = optimizer,
      innerOpt = 'trust',
      maxOuterIterations = 10,
      maxOdeRecalc = 0,
      outerMaxOdeRecalc = 0,
      print = 0,
      covMethod = '',
      calcTables = FALSE
    )
  )))
  expect_true(trialReturned)
  expect_length(values, 3)
  if (length(values) == 3) {
    expect_lt(values[3], values[2])
  }
  expect_gt(fit$env$nAnalyticGradDirect, 0)
  # The deliberately invalid trial uses the existing numerical safety path.
  expect_gt(fit$env$nFDGradFast, 0)
  expect_identical(fit$message, 'rejected trial regression')
})

test_that('FOCEI carries the likelihood score along its implicit ETA response', {
  model <- function() {
    ini({ level <- .2; error <- .2; eta ~ .3 })
    model({ d/dt(dummy) <- 0; dummy(0) <- 1
            prediction <- exp(level + eta) * dummy; prediction ~ add(error) })
  }
  data <- data.frame(ID = 1, TIME = 1:3, DV = 1.6, AMT = 0, EVID = 0)
  # The scalar score equation gives this local ETA response independently.
  pred <- exp(.2)
  mixed <- 3 * (2 * pred^2 - 1.6 * pred) / .2^2
  etaResponse <- -mixed / (mixed + 1 / .3)
  for (determinant in c('focei', 'conditional')) {
    fit <- suppressWarnings(suppressMessages(nlmixr2(
      model,
      data,
      'focei',
      foceiControl(
        fast = TRUE,
        detHessian = determinant,
        maxOuterIterations = 0,
        maxInnerIterations = 0,
        etaMat = matrix(0, 1, 1),
        resetEtaP = 0,
        print = 0,
        covMethod = '',
        calcTables = FALSE
      )
    )))
    objective <- function(level) {
      eta <- etaResponse * (level - .2)
      pred <- exp(level + eta)
      information <- if (determinant == 'focei') {
        1 / .3 + 3 * pred^2 / .2^2
      } else {
        1 / .3 + 3 * (2 * pred^2 - 1.6 * pred) / .2^2
      }
      3 * (1.6 - pred)^2 / .2^2 + eta^2 / .3 + log(information)
    }
    expected <- (objective(.20001) - objective(.19999)) / .00002
    expect_equal(fit$env$.gradDirectFirst[1], expected, tolerance = 1e-7)
    expect_gt(fit$env$nAnalyticGradDirect, 0)
    expect_equal(fit$env$nFDGradFast, 0)
  }
})
