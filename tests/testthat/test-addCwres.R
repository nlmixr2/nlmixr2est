nmTest({

  test_that("addCwres", {

    one.compartment <- function() {
      ini({
        tka <- log(1.57)
        tcl <- log(2.72)
        tv <- log(31.5)
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }

    suppressMessages(
      fitNoEta <- nlmixr2(one.compartment, theo_sd,  est="focei", control = list(print=0))
    )
    expect_true(inherits(fitNoEta$parHistData, "data.frame"))
    expect_error(
      addCwres(fitNoEta),
      regexp = "cannot add CWRES to a model without etas"
    )
  })
})
