nmTest({
  test_that("update estimates", {

    one.compartment <- function() {
      ini({
        tka <- 0.45 ; label("Log Ka")
        tcl <- 1 ; label("Log Cl")
        tv <- 3.45 ; label("Log V")
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }

    fit2 <-
      suppressMessages(suppressWarnings(
        nlmixr(
          one.compartment, theo_sd,
          est = "focei",
          control = list(print = 0)
        )
      ))

    expect_equal(fit2$theta, fit2$ui$theta)

  })
})
