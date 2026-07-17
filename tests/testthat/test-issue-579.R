nmTest({
  test_that("Issue nlmixr2est#579: saem endpoint mis-match error is actionable", {
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
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.sd)
      })
    }

    # A single-endpoint model but data with two DVID levels -> the data has more
    # observation compartments than the model has endpoints.
    d_fit <- nlmixr2data::theo_sd
    d_fit$DVID <- 1
    d_fit$DVID[d_fit$ID <= 6] <- 2

    # The error should now name the model/data endpoint counts and point the user
    # at the CMT/DVID columns instead of the old opaque message.
    expect_error(
      .nlmixr(one.compartment, d_fit, est = "saem", control = saemControlFast),
      "endpoints between the model \\(1\\) and the data \\(2\\)"
    )
    expect_error(
      .nlmixr(one.compartment, d_fit, est = "saem", control = saemControlFast),
      "CMT'/'DVID"
    )
  })
})
