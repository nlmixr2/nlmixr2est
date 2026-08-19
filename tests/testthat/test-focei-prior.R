nmTest({
  skip_if_not(exists("rxPriorBuildSpec", envir = asNamespace("rxode2"), inherits = FALSE),
              "rxode2 without the shared prior kernel (nlmixr2/rxode2#1270)")

  .oneCmt <- function() {
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
      linCmt() ~ add(add.sd)
    })
  }

  .oneCmtPrior <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
      prior(tcl) ~ dnorm(1, 0.05)
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  .oneCmtOmegaPrior <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
      prior(eta.cl) ~ dnorm(0, 0.3)
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("FOCEi's family refuses a prior that touches omega", {
    expect_error(
      suppressWarnings(suppressMessages(
        nlmixr2(.oneCmtOmegaPrior, nlmixr2data::theo_sd, est = "focei",
                control = foceiControl(maxOuterIterations = 0L, print = 0L)))),
      "omega")
  })

  test_that("the objective shifts by exactly -2*log p(theta) (#931)", {
    skip_on_cran()
    .fit0 <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, est = "posthoc")))
    .fit1 <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmtPrior, nlmixr2data::theo_sd, est = "posthoc")))
    # both fits evaluate the SAME (initial) theta, since posthoc does not
    # iterate -- so the only difference in the objective is the added prior
    # term, evaluated once at tcl's initial estimate (also the prior mean).
    .expected <- -2 * dnorm(1, 1, 0.05, log = TRUE)
    expect_equal(.fit1$objective - .fit0$objective, .expected, tolerance = 1e-6)
  })

  test_that("a prior downgrades fast=TRUE (no analytic d/dtheta log p(theta) term yet)", {
    skip_on_cran()
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmtPrior, nlmixr2data::theo_sd, est = "foceif",
              control = foceiControl(print = 0L))))
    expect_false(isTRUE(.fit$foceiControl$fast))
  })

  test_that("a prior downgrades covType='analytic' to a finite-difference covariance", {
    skip_on_cran()
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmtPrior, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(covMethod = "analytic", maxOuterIterations = 0L,
                                     print = 0L))))
    expect_false(identical(.fit$foceiControl$covType, "analytic"))
  })

  test_that("a strong prior pulls the estimate toward the prior mean (#931)", {
    skip_on_cran()
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmtPrior, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(print = 0L))))
    # tcl's prior is dnorm(1, 0.05), much tighter than the data's own
    # information about tcl -- the converged estimate should land close to
    # the prior mean, not at the (much larger) unconstrained MLE.
    expect_equal(unname(.fit$theta["tcl"]), 1, tolerance = 0.05)
  })
})
