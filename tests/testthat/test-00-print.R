test_that("print works with combined zero and correlated etas (#359)", {
  one.compartment <- function() {
    ini({
      tka <- log(1.57); label("Ka")
      tcl <- log(2.72); label("Cl")
      tv <- log(31.5); label("V")
      eta.ka ~ 0
      eta.cl + eta.v ~ c(0.3, 0.01, 0.1)
      add.sd <- 0.7
    })
    # and a model block with the error specification and model specification
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  suppressMessages(
    fit <-.nlmixr(one.compartment, theo_sd, est="saem", control = saemControlFast)
  )
  expect_output(
    .getCorPrint(fit$omegaR),
    regexp = "cor:eta.v,eta.cl",
    fixed = TRUE
  )
})
