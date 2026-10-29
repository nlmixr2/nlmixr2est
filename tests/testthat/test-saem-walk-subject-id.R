# saem's prediction walk reads a lagged dose's time through the lag of the
# subject it walks.  With fewer chain-subjects (subjects x nmc) than model
# states, the walk read the lag of a subject past the end of the solve's
# subjects, which crashed R.
.saemWalkLagModel <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    tlag <- log(0.5)
    add.sd <- 0.7
    eta.ka ~ 0.6
    eta.lag ~ 0.2
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl)
    v <- exp(tv)
    alag(depot) <- exp(tlag + eta.lag)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

nmTest({
  test_that("saem fits a lagged-dose model with fewer chain-subjects than states", {
    .d <- nlmixr2data::theo_sd[nlmixr2data::theo_sd$ID %in% 1:2, ]
    .ctl <- saemControl(nBurn = 10, nEm = 10, nmc = 1, print = 0, seed = 1L, covMethod = "", calcTables = FALSE)
    .fit <- .nlmixr(.saemWalkLagModel, .d, est = "saem", control = .ctl)
    expect_equal(.fit$objf, 76.380378, tolerance = 1e-6)
    expect_equal(
      unname(.fit$theta),
      c(1.118834, 0.731403, 3.462208, -2.023712, 0.949486),
      tolerance = 1e-5
    )
  })
})
