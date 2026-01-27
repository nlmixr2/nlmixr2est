nmTest({
  test_that("test focei", {

    one.compartment <- function() {
      ini({
        tka <- 0.45; label("Ka")
        tcl <- 1; label("Cl")
        tv <- 3.45; label("V")
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      # and a model block with the error specification and model specification
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        a <- "<5"
        if (cp >= 5) {
          a <- ">=5"
        }
        cp ~ add(add.sd)
      })
    }

    f <- .nlmixr(one.compartment, theo_sd, "focei",
                 control=foceiControl(print=0, maxOuterIterations = 1L,
                                      maxInnerIterations = 1L))

    expect_true(inherits(f$a, "factor"))

    expect_equal(unique(f$a),
                 structure(1:2, levels = c("<5", ">=5"),
                           class = "factor"))

    f <- suppressMessages(addNpde(f))

    expect_true(inherits(f$a, "factor"))

    expect_equal(unique(f$a),
                 structure(1:2, levels = c("<5", ">=5"),
                           class = "factor"))
  })
  test_that("test saem", {

    one.compartment <- function() {
      ini({
        tka <- 0.45; label("Ka")
        tcl <- 1; label("Cl")
        tv <- 3.45; label("V")
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      # and a model block with the error specification and model specification
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        a <- "<5"
        if (cp >= 5) {
          a <- ">=5"
        }
        cp ~ add(add.sd)
      })
    }



    f <- .nlmixr(one.compartment, theo_sd, "saem", control = saemControlFast)
    expect_equal(unique(f$a), factor(1:2, labels = c("<5", ">=5")))

    f <- suppressMessages(addNpde(f))
    expect_equal(unique(f$a), factor(1:2, labels = c("<5", ">=5")))

    f <- suppressMessages(addCwres(f))
    expect_equal(unique(f$a), factor(1:2, labels = c("<5", ">=5")))
  })

  test_that("test nlme", {

    one.compartment <- function() {
      ini({
        tka <- 0.45 # Log Ka
        tcl <- 1 # Log Cl
        tv <- 3.45    # Log V
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        a <- "<5"
        if (cp >= 5) {
          a <- ">=5"
        }
        cp ~ add(add.sd)
      })
    }

    f  <- .nlmixr(one.compartment, theo_sd, "nlme", control=nlmeControl(verbose=FALSE, returnObject=TRUE))

    expect_true(inherits(f$a, "factor"))

    expect_equal(unique(f$a),
                 structure(1:2, levels = c("<5", ">=5"),
                           class = "factor"))

    f <- suppressMessages(addNpde(f))

    expect_true(inherits(f$a, "factor"))

    expect_equal(unique(f$a),
                 structure(1:2, levels = c("<5", ">=5"),
                           class = "factor"))

    f <- suppressMessages(addCwres(f))

    expect_true(inherits(f$a, "factor"))

    expect_equal(unique(f$a),
                 structure(1:2, levels = c("<5", ">=5"),
                           class = "factor"))
  })

  test_that("nlm/nls", {

    one.compartment <- function() {
      ini({
        tka <- 0.45 # Log Ka
        tcl <- 1 # Log Cl
        tv <- 3.45    # Log V
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        a <- "<5"
        if (cp >= 5) {
          a <- ">=5"
        }
        cp ~ add(add.sd)
      })
    }

    f  <- .nlmixr(one.compartment, theo_sd, "nlm", nlmControl(print=0L))

    expect_true(inherits(f$a, "factor"))

    expect_equal(unique(f$a),
                 structure(1:2, levels = c("<5", ">=5"),
                           class = "factor"))

    f  <- .nlmixr(one.compartment, theo_sd, "nls", nlsControl(print=0L))

    expect_true(inherits(f$a, "factor"))

    expect_equal(unique(f$a),
                 structure(1:2, levels = c("<5", ">=5"),
                           class = "factor"))


  })


})
