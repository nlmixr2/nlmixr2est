nmTest({
  m <- function() {
    ini({
      tka <- log(1.5)
      tcl <- log(2.7)
      tv <- log(31.5)
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
      ll(err) ~ -0.5 * log(2 * pi) - log(add.sd) - 0.5 * ((DV - cp) / add.sd)^2
    })
  }
  # The Lagrangian's eta-stationarity term: the inner mode is solved to a tolerance
  # (n1qn1 epsilon), so the gradient must not carry it.  Without the term the mu-theta
  # rows moved ~2e-3 between the default and a converged inner mode.
  test_that("forward ll() fast gradient does not carry the inner tolerance", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .fit <- function(eps) {
      suppressMessages(nlmixr2(
        m,
        nlmixr2data::theo_sd,
        "focei",
        foceiControl(
          print = 0L,
          covMethod = "",
          fast = TRUE,
          maxOuterIterations = 0L,
          epsilon = eps,
          rxControl = rxode2::rxControl(
            atol = 1e-11,
            rtol = 1e-11,
            atolSens = 1e-11,
            rtolSens = 1e-11,
            method = "dop853"
          )
        )
      ))
    }
    fL <- .fit(1e-3)
    fT <- .fit(1e-12)
    th <- fL$theta * 1.02
    om <- fL$omega * 1.1
    gL <- nlmixr2ObjGrad(fL, theta = th, omega = om)
    gT <- nlmixr2ObjGrad(fT, theta = th, omega = om)
    expect_gt(attr(gL, "nAnalyticGradDirect"), 0L)
    expect_equal(gL$gradient, gT$gradient, tolerance = 2e-5)
  })
})
