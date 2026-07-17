# SAEM covMethod="analytic": the FOCEI analytic observed-information covariance
# at the converged SAEM estimates, falling back to the linearized FIM (linFim)
# when out of analytic scope.  Weekly batch (multi-iteration fits).

nmTest({
  odeMod <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  linMod <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  ctl <- function(...) saemControl(nBurn = 100, nEm = 150, print = 0, seed = 42, ...)

  test_that("SAEM covMethod='analytic' installs the analytic covariance", {
    skip_on_cran()
    ## the default is "sa"; request the analytic observed information explicitly
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(odeMod, nlmixr2data::theo_sd, est = "saem",
              control = ctl(covMethod = "analytic"))))
    expect_identical(fit$covMethod, "analytic")
    expect_true(all(is.finite(fit$parFixedDf$SE)))
    expect_true(all(fit$parFixedDf$SE > 0))
    ## the linFim fallback is retained and selectable
    expect_true("linFim" %in% names(fit$env$covList))
    expect_error(setCov(fit, "linFim"), NA)
    expect_identical(fit$covMethod, "linFim")
  })

  test_that("explicit covMethod='linFim' skips the analytic attempt", {
    skip_on_cran()
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(odeMod, nlmixr2data::theo_sd, est = "saem",
              control = ctl(covMethod = "linFim"))))
    expect_identical(fit$covMethod, "linFim")
  })

  test_that("out-of-scope model (linCmt) with covMethod='analytic' falls back to linFim", {
    skip_on_cran()
    expect_message(
      fit <<- suppressWarnings(
        nlmixr2(linMod, nlmixr2data::theo_sd, est = "saem",
                control = ctl(covMethod = "analytic"))),
      "linearized FIM")
    expect_identical(fit$covMethod, "linFim")
    expect_true(all(is.finite(fit$parFixedDf$SE)))
  })

  test_that("saem default covMethod is the stochastic-approximation FIM", {
    skip_on_cran()
    expect_identical(saemControl()$covMethod, "sa")
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(odeMod, nlmixr2data::theo_sd, est = "saem", control = ctl())))
    expect_identical(fit$covMethod, "sa")
    expect_true(all(is.finite(fit$parFixedDf$SE)))
  })

  test_that("fsaem produces a labeled covariance with the sa default", {
    skip_on_cran()
    ## fsaem is saemControl(fast=TRUE): it inherits the "sa" covMethod default.
    ## The fast kernel overwrites the ui covMethod when it sets up its FOCEi inner
    ## problem; .saemFamilyFit restores it so the covariance is computed + labeled.
    expect_identical(fsaemControl()$covMethod, "sa")
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(odeMod, nlmixr2data::theo_sd, est = "fsaem", control = ctl())))
    expect_s3_class(fit, "nlmixr2FitData")
    expect_true(is.finite(fit$objf))
    expect_identical(fit$covMethod, "sa")
    expect_false(is.null(fit$cov))
    expect_true(all(is.finite(fit$parFixedDf$SE)))
    expect_true(all(fit$parFixedDf$SE > 0))
  })
})
