nmTest({
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

  fit <- .nlmixr(one.compartment, theo_sd, est = "saem", control = saemControlFast)

  test_that("as.rxUi works for estimated models", {
    expect_s3_class(as.rxUi(fit), "rxUi")
  })
})
