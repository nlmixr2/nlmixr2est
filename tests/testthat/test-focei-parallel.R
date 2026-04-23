nmTest({
  test_that("FOCEi with cores=2 matches cores=1 result", {
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
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv  + eta.v)
        d/dt(depot)  <- -ka * depot
        d/dt(center) <-  ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    fit1 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "",
                                          rxControl = rxode2::rxControl(cores = 1L)))
    fit2 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "",
                                          rxControl = rxode2::rxControl(cores = 2L)))
    expect_s3_class(fit1, "nlmixr2FitData")
    expect_s3_class(fit2, "nlmixr2FitData")
    expect_equal(fit1$objective, fit2$objective, tolerance = 1e-4)
    expect_equal(as.numeric(fit1$theta), as.numeric(fit2$theta), tolerance = 1e-4)
  })
})
