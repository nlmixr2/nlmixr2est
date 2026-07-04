nmTest({
  test_that("single population parameter estimation becomes optimize", {

    one.compartment <- function() {
      ini({
        tka <- fix(0.45)
        tcl <- 1
        tv <- fix(3.45)
        eta.ka ~ fix(0.6)
        eta.cl ~ fix(0.3)
        eta.v ~ fix(0.1)
        add.sd <- fix(0.7)
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.sd)
      })
    }

    fit <-
      .nlmixr(
        one.compartment, theo_sd,
        est="focei",
        control = foceiControl(print=0)
      )
      expect_equal(fit$message, "stats::optimize for 1 dimensional optimization")
  })

  test_that("focei mixed effects works on matExp", {
    one.compartment <- function() {
      ini({
        tka <- 0.45
        tcl <- log(c(0, 2.7, 100))
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        matExp()
        k_depot_central = exp(tka + eta.ka)
        k_central_output = exp(tcl + eta.cl) / exp(tv + eta.v)
        cp = central / exp(tv + eta.v)
        cp ~ add(add.sd)
      })
    }

    s <- rxUiGet.foceiEtaS(list(rxode2::rxode2(one.compartment)))
    expect_true(exists("..jacobian", envir = s, inherits = FALSE))
    expect_true(exists("..sens", envir = s, inherits = FALSE))
  })

  test_that("focei matExp event sensitivities build HdEta", {
    one.compartment <- function() {
      ini({
        tka <- 0.45
        tf <- 0
        eta.ka ~ 0.3
        eta.f ~ 0.2
        add.sd <- 0.7
      })
      model({
        matExp()
        k_depot_central = exp(tka + eta.ka)
        k_central_output = 0.2
        f(depot) <- expit(tf + eta.f)
        cp = central / 10
        cp ~ add(add.sd)
      })
    }

    ui <- rxode2::rxode2(one.compartment)

    s <- rxUiGet.foceiHdEta(list(ui))
    expect_true(exists("..HdEta", envir = s, inherits = FALSE))
    expect_true(any(grepl("rx__sens_central_BY_ETA_2__", get("..HdEta", envir = s))))
  })

  test_that("all parameters fixed becomes posthoc", {

    one.compartment <- function() {
      ini({
        tka <- fix(0.45)
        tcl <- fix(1)
        tv <- fix(3.45)
        eta.ka ~ fix(0.6)
        eta.cl ~ fix(0.3)
        eta.v ~ fix(0.1)
        add.sd <- fix(0.7)
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.sd)
      })
    }

    fit <- .nlmixr(one.compartment, theo_sd, est="focei", control=list(print=0))
    expect_s3_class(fit, "nlmixr2FitData")
  })

  test_that("single population parameter estimation becomes optimize (no eta)", {
    one.compartment <- function() {
      ini({
        tka <- fix(0.45)
        tcl <- 1
        tv <- fix(3.45)
        add.sd <- fix(0.7)
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.sd)
      })
    }

    fit <- .nlmixr(one.compartment, theo_sd, est = "focei", control=list(print=0))
    expect_equal(fit$message, "stats::optimize for 1 dimensional optimization")
  })

  test_that("all parameters fixed with no etas errors", {
    # nothing to estimate, should error
    one.compartment <- function() {
      ini({
        tka <- fix(0.45)
        tcl <- fix(1)
        tv <- fix(3.45)
        add.sd <- fix(0.7)
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.sd)
      })
    }
    expect_error(
      .nlmixr(one.compartment, theo_sd, est = "focei", control=list(print=0)),
      regexp = "no parameters to estimate"
    )
  })
})
