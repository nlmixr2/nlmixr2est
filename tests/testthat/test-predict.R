nmTest({
  test_that("predict", {

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
        cp ~ add(add.sd)
      })
    }

    ## The fit is performed by the function nlmixr/nlmix2 specifying the model, data and estimate
    fit <- nlmixr2(one.compartment, theo_sd,  est="focei",
                   foceiControl(maxOuterIterations = 0L))

    md <- do.call("predict", c(list(fit, theo_md), fit$control))

    md2 <- nlmixr2(fit, theo_md, "predict")

    expect_equal(as.data.frame(md), as.data.frame(md2), tolerance = 1e-4)

    md <- predict(fit, theo_md)

    expect_equal(as.data.frame(md), as.data.frame(md2), tolerance = 1e-4)

  })
})
