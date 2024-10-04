nmTest({
  .nlmixr <- function(...) suppressWarnings(suppressMessages(nlmixr(...)))

  test_that("single population parameter estimation becomes optimize", {

    one.compartment <- function() {
      ini({
        tka <- fix(0.45) # Log Ka
        tcl <- 1 # Log Cl
        tv <- fix(3.45)    # Log V
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

    fit <- .nlmixr(one.compartment, theo_sd,  est="focei", control=list(print=0))
    expect_equal(fit$message, "stats::optimize for 1 dimensional optimization")
  })
  test_that("all parameters fixed becomes posthoc", {

    one.compartment <- function() {
      ini({
        tka <- fix(0.45) # Log Ka
        tcl <- fix(1) # Log Cl
        tv <- fix(3.45)    # Log V
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

    fit <- .nlmixr(one.compartment, theo_sd,  est="focei", control=list(print=0))
    expect_true(inherits(fit, "nlmixr2FitData"))
  })

  test_that("single population parameter estimation becomes optimize (no eta)", {

    one.compartment <- function() {
      ini({
        tka <- fix(0.45) # Log Ka
        tcl <- 1 # Log Cl
        tv <- fix(3.45)    # Log V
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

    fit <- .nlmixr(one.compartment, theo_sd,  est="focei", control=list(print=0))
    expect_equal(fit$message, "stats::optimize for 1 dimensional optimization")
  })

  test_that("all parameters fixed with no etas errors", {
    # nothing to estimate, should error
    one.compartment <- function() {
      ini({
        tka <- fix(0.45) # Log Ka
        tcl <- fix(1) # Log Cl
        tv <- fix(3.45)    # Log V
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
    expect_error(.nlmixr(one.compartment, theo_sd,  est="focei", control=list(print=0)))
  })
})
