# An M3/M4 row carries no log(2*pi) on either likelihood route, as in NONMEM (#1478)
test_that("censored OFV constant matches NONMEM on both likelihood routes (#1478)", {
  skip_on_cran()
  set.seed(42)
  tt <- c(0.5, 1, 2, 4, 8)
  d <- do.call(
    rbind,
    lapply(1:40, function(i) {
      a <- 10 * exp(rnorm(1, 0, sqrt(0.1)))
      k <- 0.3 * exp(rnorm(1, 0, sqrt(0.1)))
      data.frame(ID = i, TIME = tt, DV = a * exp(-k * tt) + 0.5 * rnorm(5))
    })
  )
  d$CENS <- as.integer(d$DV < 3)
  d$DV[d$CENS == 1] <- 3
  nUncens <- sum(d$CENS == 0)
  expect_equal(nUncens, 143L)

  gauss <- function() {
    ini({
      ta <- log(9)
      tk <- log(0.25)
      ab <- 0.7
      eta.a ~ 0.15
      eta.k ~ 0.2
    })
    model({
      a <- exp(ta + eta.a)
      k <- exp(tk + eta.k)
      f <- a * exp(-k * time)
      f ~ add(ab)
    })
  }
  llik <- function() {
    ini({
      ta <- log(9)
      tk <- log(0.25)
      ab <- 0.7
      eta.a ~ 0.15
      eta.k ~ 0.2
    })
    model({
      a <- exp(ta + eta.a)
      k <- exp(tk + eta.k)
      f <- a * exp(-k * time)
      f ~ add(ab) + dnorm()
    })
  }
  .ctl <- list(maxOuterIterations = 0, covMethod = "", print = 0, calcTables = FALSE)
  fitGauss <- .nlmixr(gauss, d, "flaplace", control = .ctl)
  fitLlik <- .nlmixr(llik, d, "laplace", control = .ctl)
  # NONMEM 7.6.0 METHOD=COND LAPLACE INTER MAXEVAL=0 on the same data: 213.55289
  expect_equal(fitGauss$objf, 213.55289, tolerance = 1e-5)
  expect_equal(fitLlik$objf, 213.55289, tolerance = 1e-5)
  # -2LL adds log(2*pi) for the uncensored rows only (NONMEM LNTWOPI: 476.36931)
  for (.fit in list(fitGauss, fitLlik)) {
    expect_equal(-2 * as.numeric(logLik(.fit)) - .fit$objf, nUncens * log(2 * pi), tolerance = 1e-8)
    expect_equal(-2 * as.numeric(logLik(.fit)), 476.36931, tolerance = 1e-5)
  }
})
