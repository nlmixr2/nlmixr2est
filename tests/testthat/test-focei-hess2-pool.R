test_that('conditional inner Hessians use their own event shape in an outer-sized pool', {
  model <- function() {
    ini({ level <- log(1.5); logk <- log(.7); error <- .2; eta ~ .3 })
    model({ k <- exp(logk); central(0) <- exp(level + eta)
            d/dt(central) <- -k * central; prediction <- central
            prediction ~ add(error) })
  }
  times <- c(.5, 1, 2, 4)
  base <- 1.5 * exp(-.7 * times)
  data <- data.frame(ID = rep(1:2, each = 4), TIME = rep(times, 2), DV = c(1.1 * base, .9 * base), AMT = 0, EVID = 0)
  expected <- vapply(
    1:2,
    function(id) {
      y <- data$DV[data$ID == id]
      uniroot(
        function(eta) {
          pred <- base * exp(eta)
          sum((pred - y) * pred) / .2^2 + eta / .3
        },
        c(-1, 1),
        tol = 1e-12
      )$root
    },
    numeric(1)
  )
  fit <- suppressWarnings(suppressMessages(nlmixr2(
    model,
    data,
    'focei',
    foceiControl(
      fast = TRUE,
      innerOpt = 'trust',
      innerHessian = 'conditional',
      maxOuterIterations = 0,
      maxInnerIterations = 100,
      sigdig = 6,
      trustFterm = 1e-12,
      trustMterm = 1e-12,
      etaNudge = 0,
      etaNudge2 = 0,
      rxControl = rxode2::rxControl(maxsteps = 2000L, atol = 1e-10, rtol = 1e-8),
      print = 0,
      covMethod = '',
      calcTables = FALSE
    )
  )))
  info <- .odeSwapInfo()
  hess2 <- info$models$neq[match('hess2', info$models$name)]
  expect_gt(info$poolNeq, hess2)
  expect_equal(info$overrideNeutralizedN, 0)
  expect_equal(as.numeric(fit$eta[[2]]), expected, tolerance = 1e-5)
  expect_gt(fit$env$nConditionalInnerHessian, 0)
  expect_gt(fit$env$nAnalyticGradDirect, 0)
  expect_equal(fit$env$nFDGradFast, 0)
})


test_that('conditional Hessians preserve all mixed derivatives with three ETAs', {
  model <- function() {
    ini({ logk <- log(.7); error <- .2; eta1 ~ .3; eta2 ~ .2; eta3 ~ .1 })
    model({ k <- exp(logk); central(0) <- 1.5 * exp(eta1 + 2 * eta2 + 3 * eta3)
            d/dt(central) <- -k * central; prediction <- central
            prediction ~ add(error) })
  }
  times <- c(.5, 1, 2, 4)
  base <- 1.5 * exp(-.7 * times)
  data <- data.frame(ID = rep(1:2, each = 4), TIME = rep(times, 2), DV = c(1.1 * base, .9 * base), AMT = 0, EVID = 0)
  variance <- c(.3, .2, .1)
  weights <- 1:3
  totalVariance <- sum(variance * weights^2)
  expected <- t(vapply(
    1:2,
    function(id) {
      y <- data$DV[data$ID == id]
      z <- uniroot(
        function(z) {
          pred <- base * exp(z)
          sum((pred - y) * pred) / .2^2 + z / totalVariance
        },
        c(-1, 1),
        tol = 1e-12
      )$root
      variance * weights / totalVariance * z
    },
    numeric(3)
  ))
  fit <- suppressWarnings(suppressMessages(nlmixr2(
    model,
    data,
    'focei',
    foceiControl(
      fast = TRUE,
      innerOpt = 'trust',
      innerHessian = 'conditional',
      maxOuterIterations = 0,
      maxInnerIterations = 100,
      sigdig = 6,
      trustFterm = 1e-12,
      trustMterm = 1e-12,
      etaNudge = 0,
      etaNudge2 = 0,
      rxControl = rxode2::rxControl(maxsteps = 2000L, atol = 1e-10, rtol = 1e-8),
      print = 0,
      covMethod = '',
      calcTables = FALSE
    )
  )))
  expect_equal(unname(as.matrix(fit$eta[-1])), expected, tolerance = 1e-5)
  expect_equal(.odeSwapInfo()$overrideNeutralizedN, 0)
  expect_gt(fit$env$nAnalyticGradDirect, 0)
  expect_equal(fit$env$nFDGradFast, 0)
})

test_that('generalized inner Hessians share the correct solver role', {
  model <- function() {
    ini({ level <- log(4); logk <- log(.7); eta ~ .3 })
    model({ k <- exp(logk); central(0) <- exp(level + eta)
            d/dt(central) <- -k * central; effect <- central
            effect ~ dpois(effect) })
  }
  times <- c(.5, 1, 2, 4)
  base <- 4 * exp(-.7 * times)
  data <- data.frame(ID = rep(1:2, each = 4), TIME = rep(times, 2), DV = c(3, 2, 1, 0, 1, 1, 1, 0), AMT = 0, EVID = 0)
  expected <- vapply(
    1:2,
    function(id) {
      y <- data$DV[data$ID == id]
      uniroot(function(eta) sum(base * exp(eta) - y) + eta / .3, c(-2, 2), tol = 1e-12)$root
    },
    numeric(1)
  )
  fit <- suppressWarnings(suppressMessages(nlmixr2(
    model,
    data,
    'focei',
    foceiControl(
      fast = TRUE,
      innerOpt = 'n1qn1',
      warm = 'calc',
      epsilon = 1e-8,
      maxOuterIterations = 0,
      maxInnerIterations = 100,
      sigdig = 6,
      etaNudge = 0,
      etaNudge2 = 0,
      rxControl = rxode2::rxControl(maxsteps = 2000L, atol = 1e-10, rtol = 1e-8),
      print = 0,
      covMethod = '',
      calcTables = FALSE
    )
  )))
  expect_equal(as.numeric(fit$eta[[2]]), expected, tolerance = 1e-5)
  expect_equal(.odeSwapInfo()$overrideNeutralizedN, 0)
  expect_gt(fit$env$nAnalyticGradDirect, 0)
  expect_equal(fit$env$nFDGradFast, 0)
})
