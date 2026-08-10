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

  test_that("the FIM's single residual slot is filled only when it fits", {
    skip_on_cran(); skip_if_not_installed("nlmixr2data")
    # The FIM carries exactly ONE residual slot (nb_param = nphi1 + nlambda + 1),
    # and d1_logsigma2 = 0.5*resy/sigma2 - 0.5*ntotal counts EVERY observation
    # against endpoint 0's variance.  That only describes a model with one
    # estimated residual, additive.  Anything else used to fill the slot anyway
    # and invert it together with the parameters that ARE reported.
    #
    # Two normal endpoints, residual SDs an order of magnitude apart -- the case
    # where mixing one endpoint's variance with another's residual is worst.
    .d <- mkPkpdTwoEndpointData(n = 20L)
    .slotZero <- function(f) {
      .H <- f$env$saem$Ha
      .n <- nrow(.H)
      all(.H[.n, ] == 0) && all(.H[, .n] == 0)
    }
    .fit <- function(m, dd) suppressMessages(suppressWarnings(
      nlmixr2(m, dd, est = "saem",
              control = saemControl(nBurn = 60, nEm = 30, nmc = 3, seed = 7,
                                    print = 0L, covMethod = "fim"))))

    .two <- .fit(pkpd.two.endpoint, .d)
    expect_true(.slotZero(.two))
    # and the method still stands up: dropping the row happens BEFORE the
    # inversion, so theta + Omega come out rather than the whole FIM being lost
    # to a singular matrix
    expect_equal(.two$covMethod, "fim")
    expect_true(all(is.finite(stats::na.omit(.two$parFixedDf$SE))))

    # a single additive endpoint is the one case the slot describes: still filled
    .pk <- .d[.d$dvid == "cp" | .d$evid == 1, ]
    .one <- function() {
      ini({
        tcl <- -3.2; tv <- -1
        eta.cl ~ 0.09; eta.v ~ 0.09
        pk.sd <- 0.4
      })
      model({
        ka <- exp(0.5); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(pk.sd)
      })
    }
    expect_false(.slotZero(.fit(.one, .pk)))

    # a proportional residual is not what the slot describes either: .saemFimToCov
    # never reported it, but it used to sit in the matrix being inverted
    .prop <- function() {
      ini({
        tcl <- -3.2; tv <- -1
        eta.cl ~ 0.09; eta.v ~ 0.09
        pk.p <- 0.1
      })
      model({
        ka <- exp(0.5); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ prop(pk.p)
      })
    }
    .pr <- .fit(.prop, .pk)
    expect_true(.slotZero(.pr))
    expect_equal(.pr$covMethod, "fim")
  })

  test_that("covMethod='fim' is refused for a general log-likelihood endpoint", {
    skip_on_cran(); skip_if_not_installed("nlmixr2data")
    # An ll() endpoint has no residual parameter, so the slot is already zeroed and
    # dropping it before inverting would make the FIM invertible.  Measured against
    # an exact marginal likelihood that answer is several-fold too small (#871), so
    # the combination is refused up front and the linearized FIM is used instead.
    .d <- mkPkpdTwoEndpointData(n = 20L)
    .pk <- .d[.d$dvid == "cp" | .d$evid == 1, ]
    .llm <- function() {
      ini({
        tcl <- -3.2; tv <- -1
        eta.cl ~ 0.09; eta.v ~ 0.09
        lsd <- log(0.4)
      })
      model({
        ka <- exp(0.5); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        s <- exp(lsd)
        ll(cp) ~ -0.5 * log(2 * pi) - lsd - 0.5 * ((DV - cp) / s)^2
      })
    }
    .f <- suppressMessages(suppressWarnings(
      nlmixr2(.llm, .pk, est = "saem",
              control = saemControl(nBurn = 60, nEm = 30, nmc = 3, seed = 7,
                                    print = 0L, covMethod = "fim"))))
    expect_equal(.f$covMethod, "linFim")
    expect_true(all(is.finite(stats::na.omit(.f$parFixedDf$SE))))
  })
})
