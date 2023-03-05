nmTest({
  test_that("boundary issue", {

    one.compartment <- function() {
      ini({
        tka <- c(-6, -4, 2)
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)*100
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.sd)
      })
    }

    ## The fit is performed by the function nlmixr/nlmix2 specifying the model, data and estimate
    fit <- nlmixr2(one.compartment, theo_sd,  est="focei", control = list(print=0))
    expect_true(inherits(fit$cov, "matrix"))

  })
})
