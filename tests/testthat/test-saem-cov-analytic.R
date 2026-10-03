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
      nlmixr2(odeMod, nlmixr2data::theo_sd, est = "saem", control = ctl(covMethod = "analytic"))
    ))
    expect_identical(.covBaseName(fit$covMethod), "analytic")
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
      nlmixr2(odeMod, nlmixr2data::theo_sd, est = "saem", control = ctl(covMethod = "linFim"))
    ))
    expect_identical(fit$covMethod, "linFim")
  })

  test_that("out-of-scope model (linCmt) with covMethod='analytic' falls back to linFim", {
    skip_on_cran()
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(linMod, nlmixr2data::theo_sd, est = "saem", control = ctl(covMethod = "analytic"))
    ))
    # a warning, so it is kept in $runInfo
    expect_true("\"analytic\" covariance could not be computed; kept \"linFim\"" %in% fit$runInfo)
    expect_identical(fit$covMethod, "linFim")
    expect_true(all(is.finite(fit$parFixedDf$SE)))
  })

  test_that("SAEM covMethod='analytic' uses the FOCEI formulas of SAEM's likelihood", {
    skip_on_cran()
    # with a proportional error the interaction changes the observed information; the
    # output step that finalizes a saem fit evaluates FOCE at SAEM's ETAs, and the
    # analytic covariance used to take that interaction = 0 from fit$finalUi
    ceMod <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.3; prop.sd <- 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd) + prop(prop.sd)
      })
    }
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(ceMod, nlmixr2data::theo_sd, est = "saem", control = ctl(covMethod = "analytic"))
    ))
    expect_identical(.covBaseName(fit$covMethod), "analytic")
    expect_identical(rxode2::rxGetControl(fit$finalUi, "interaction", NA), 0L)
    expect_identical(.foceiFitInteraction(fit, fit$finalUi), 1L)
    # the analytic covariance a FOCEI (and a FOCE) fit assembles at SAEM's estimates and ETAs
    .ref <- function(interaction) {
      suppressMessages(suppressWarnings(nlmixr2(
        fit$finalUi,
        nlmixr2data::theo_sd,
        est = "focei",
        control = foceiControl(
          print = 0,
          maxOuterIterations = 0L,
          maxInnerIterations = 0L,
          etaMat = fit$etaMat,
          covMethod = "analytic",
          interaction = interaction
        )
      )))
    }
    .se <- function(f) sqrt(diag(f$cov))
    .focei <- .ref(TRUE)
    .foce <- .ref(FALSE)
    expect_equal(.se(fit), .se(.focei), tolerance = 1e-6)
    # the add.sd SE is 45% larger with the FOCE formulas here
    expect_gt(.se(.foce)[["add.sd"]] / .se(.focei)[["add.sd"]], 1.2)
    # setCov(fit, "analytic") assembles the same on a fit that has no analytic covariance yet
    fit2 <- suppressMessages(suppressWarnings(
      nlmixr2(ceMod, nlmixr2data::theo_sd, est = "saem", control = ctl(covMethod = "linFim"))
    ))
    suppressMessages(suppressWarnings(setCov(fit2, "analytic")))
    expect_identical(fit2$covMethod, "analytic")
    expect_equal(.se(fit2), .se(.focei)[names(.se(fit2))], tolerance = 1e-6)
  })

  test_that("saem default covMethod is the stochastic-approximation FIM", {
    skip_on_cran()
    expect_identical(saemControl()$covMethod, "sa")
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(odeMod, nlmixr2data::theo_sd, est = "saem", control = ctl())
    ))
    expect_identical(fit$covMethod, "sa")
    expect_true(all(is.finite(fit$parFixedDf$SE)))
  })
})
