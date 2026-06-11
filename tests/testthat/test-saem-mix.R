nmTest({
  test_that("test SAEM mixture model estimation", {
    set.seed(42)
    nSubj <- 30
    subPop <- rbinom(nSubj, 1, 0.6) + 1 # 1 or 2
    clSim <- ifelse(subPop == 1, 1.2, 6.0)

    simData <- do.call(rbind, lapply(1:nSubj, function(i) {
      subjCl <- clSim[i]
      times <- c(0.5, 1, 2, 4, 8, 12, 24)
      kaVal <- 1.5
      vVal <- 24.0
      kVal <- subjCl / vVal
      cp <- 100 * kaVal / (vVal * (kaVal - kVal)) * (exp(-kVal * times) - exp(-kaVal * times)) + rnorm(length(times), 0, 0.05)
      if (any(cp < 0)) cp[cp < 0] <- 0
      data.frame(
        ID = i,
        TIME = c(0, times),
        AMT = c(100, rep(0, length(times))),
        EVID = c(1, rep(0, length(times))),
        DV = c(0, cp),
        CMT = c(1, rep(2, length(times)))
      )
    }))

    oneCompartmentMix <- function() {
      ini({
        tka <- log(1.5)
        tcl1 <- log(1.0)
        tcl2 <- log(5.0)
        tv <- log(20)
        p1 <- 0.5
        eta.cl ~ 0.01
        eta.v ~ 0.01
        eta.ka ~ 0.01
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }

    # Run SAEM estimation with mixture
    fitSaem <- expect_error(
      .nlmixr(oneCompartmentMix, simData, est="saem",
              saemControl(print = 0, seed = 1234, nBurn = 5, nEm = 5,
                          calcTables = FALSE, covMethod = 0L)),
      NA
    )
    
    # Check that estimated mixture parameter p1 exists
    expect_true("p1" %in% names(fitSaem$theta))
    
    # Check that estimated parameters are accessible
    expect_true(is.numeric(fitSaem$theta["p1"]))
  })
})
