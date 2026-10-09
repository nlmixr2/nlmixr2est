nmTest({
  # A lagged first dose sorted after a time-zero observation, which then read
  # the dose's CMT and returned the wrong endpoint's IPRED (#1182)
  test_that("time-zero IPRED uses its own endpoint with a lagged dose", {
    d <- data.frame(
      ID = 1,
      TIME = c(0, 0, 0.5, 2),
      AMT = c(100, 0, 0, 0),
      DV = c(0, 100, 100, 6),
      DVID = c(0, 2, 2, 1),
      EVID = c(1, 0, 0, 0),
      MDV = c(1, 0, 0, 0)
    )
    f <- function() {
      ini({
        tlag <- fixed(1)
        ka <- fixed(1)
        vc <- fixed(10)
        cl <- fixed(0.1)
        r0 <- fixed(100)
        kout <- fixed(0.1)
        sigmaPk <- 1
        sigmaPd <- 1
      })
      model({
        alag(depot) <- tlag
        response(0) <- r0
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl * central / vc
        d/dt(response) <- kout * r0 - kout * response
        yPk <- central / vc
        yPk ~ add(sigmaPk)
        yPd <- response
        yPd ~ add(sigmaPd)
      })
    }
    fit <- .nlmixr(f, d, est = "focei", control = foceiControl(print = 0))
    expect_equal(fit$TIME, c(0, 0.5, 2))
    expect_equal(fit$IPRED[1:2], c(100, 100))
    expect_equal(fit$IPRED[1:2], fit$yPd[1:2])
  })
})
