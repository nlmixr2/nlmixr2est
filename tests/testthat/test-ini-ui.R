nmTest({
  test_that("ini ui works", {

    one.compartment <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      # model block
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

    f <-  one.compartment()

    ## fit the model
    fit <- nlmixr2(one.compartment, theo_sd,  est="saem", saemControl(print=0))

    expect_equal(fit$iniUi$iniDf, f$iniDf)

  })
})
