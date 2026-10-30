test_that('ETA extrapolation uses the gradient-point baseline across line-search trials', {
  model <- function() {
    ini({ theta <- .2; slope <- .1; error <- .2; eta ~ .3 })
    model({ d/dt(dummy) <- 0; dummy(0) <- 1
            prediction <- exp(theta + eta) * dummy + slope * TIME; prediction ~ add(error) })
  }
  data <- data.frame(ID=rep(1:2, each=3), TIME=rep(1:3, 2),
                     DV=rep(exp(.2) + .1 * (1:3), 2), AMT=0, EVID=0)
  values <- NULL
  optimizer <- function(par, fn, gr, ...) {
    fn(par)
    gr(par)
    a <- b <- par
    a[1] <- a[1] + .8
    b[1] <- b[1] - .6
    first <- fn(a)
    fn(b)
    values <<- c(first, fn(a))
    list(x=par, convergence=0L, message='line-search regression')
  }
  # One inner step exposes the starting point; fully converged solves can hide it.
  fit <- suppressWarnings(suppressMessages(nlmixr2(
    model, data, 'focei', foceiControl(
      fast=TRUE, outerOpt=optimizer, innerOpt='trust', maxOuterIterations=10,
      maxInnerIterations=1, mceta=-2, resetEtaP=0, etaNudge=0, etaNudge2=0,
      print=0, covMethod='', calcTables=FALSE
    )
  )))
  expect_length(values, 2)
  expect_equal(values[1], values[2], tolerance=1e-10)
  expect_gt(fit$env$nAnalyticGradDirect, 0)
  expect_equal(fit$env$nFDGradFast, 0)
})
