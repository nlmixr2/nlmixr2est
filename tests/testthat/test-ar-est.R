nmTest({
  test_that("ar(1) autocorrelated residual estimation recovers cor (saem + focei)", {
    skip_on_cran()

    # simulate AR(1) residual data (continuous-time, phi = cor^dt) at cor = 0.7
    .m <- function() {
      ini({
        tka <- log(1.2)
        tcl <- log(0.2)
        tv <- log(5)
        eta.cl ~ 0.09
        add.sd <- 0.4
        ar1.cor <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp ~ add(add.sd) + ar(ar1.cor)
      })
    }
    .ev <- rxode2::et(rxode2::et(amt = 100), seq(0.5, 24, by = 1))
    .sim <- rxode2::rxWithSeed(2026, rxode2::rxSolve(.m, .ev, nSub = 40,
                                     returnType = "data.frame", addDosing = TRUE))
    .idc <- if ("sim.id" %in% names(.sim)) "sim.id" else "id"
    .dat <- data.frame(ID = .sim[[.idc]], TIME = .sim$time,
                       EVID = ifelse(is.na(.sim$evid), 0, .sim$evid),
                       AMT = ifelse(is.na(.sim$amt), 0, .sim$amt),
                       DV = ifelse(is.na(.sim$evid) | .sim$evid == 0, .sim$sim, NA))
    .dat$EVID[.dat$AMT > 0] <- 1

    # fit with cor started well away from the truth (0.3) to exercise recovery
    .mF <- function() {
      ini({
        tka <- log(1)
        tcl <- log(0.3)
        tv <- log(7)
        eta.cl ~ 0.1
        add.sd <- 0.6
        ar1.cor <- 0.3
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp ~ add(add.sd) + ar(ar1.cor)
      })
    }

    # SAEM: whitened-conditional E-step/M-step (grid-search cor update)
    .fs <- .nlmixr(.mF, .dat, est = "saem",
                   control = saemControl(nBurn = 150, nEm = 150, print = 0, seed = 42))
    .corS <- unname(.fs$parFixedDf["ar1.cor", "Estimate"])
    .sdS <- unname(.fs$parFixedDf["add.sd", "Estimate"])
    expect_equal(.corS, 0.7, tolerance = 0.12)
    expect_equal(.sdS, 0.4, tolerance = 0.12)

    # FOCEI: norm-form (mean/variance) whitened residual -> exact eta-Hessian
    .ff <- .nlmixr(.mF, .dat, est = "focei", control = foceiControl(print = 0))
    .corF <- unname(.ff$parFixedDf["ar1.cor", "Estimate"])
    .sdF <- unname(.ff$parFixedDf["add.sd", "Estimate"])
    expect_equal(.corF, 0.7, tolerance = 0.15)
    expect_equal(.sdF, 0.4, tolerance = 0.15)
    # regression guard: cor must not collapse to the [0,1) boundary (nor the sd)
    expect_lt(.corF, 0.95)
    expect_gt(.sdF, 0.2)
  })

  test_that("ar() combines with the log-likelihood residuals (+dnorm/+dt/+dcauchy) under focei", {
    skip_on_cran()

    # smaller/sparser data -- these paths just need to build and run (they
    # degenerate to the explicit log-likelihood inner model, not the norm form)
    .m <- function() {
      ini({
        tka <- log(1.2); tcl <- log(0.2); tv <- log(5)
        eta.cl ~ 0.09; add.sd <- 0.4; ar1.cor <- 0.7
      })
      model({
        ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp ~ add(add.sd) + ar(ar1.cor)
      })
    }
    .ev <- rxode2::et(rxode2::et(amt = 100), seq(0.5, 24, by = 1))
    .sim <- rxode2::rxWithSeed(2026, rxode2::rxSolve(.m, .ev, nSub = 40,
                                     returnType = "data.frame", addDosing = TRUE))
    .idc <- if ("sim.id" %in% names(.sim)) "sim.id" else "id"
    .dat <- data.frame(ID = .sim[[.idc]], TIME = .sim$time,
                       EVID = ifelse(is.na(.sim$evid), 0, .sim$evid),
                       AMT = ifelse(is.na(.sim$amt), 0, .sim$amt),
                       DV = ifelse(is.na(.sim$evid) | .sim$evid == 0, .sim$sim, NA))
    .dat$EVID[.dat$AMT > 0] <- 1

    # +dnorm(): normal via the log-likelihood path -- must match the norm form
    .mN <- function() {
      ini({
        tka <- log(1); tcl <- log(0.3); tv <- log(7)
        eta.cl ~ 0.1; add.sd <- 0.6; ar1.cor <- 0.3
      })
      model({
        ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp ~ add(add.sd) + ar(ar1.cor) + dnorm()
      })
    }
    .fN <- .nlmixr(.mN, .dat, est = "focei", control = foceiControl(print = 0))
    expect_equal(unname(.fN$parFixedDf["ar1.cor", "Estimate"]), 0.7, tolerance = 0.15)
    expect_true(is.finite(.fN$objf))

    # +dt(df) and +dcauchy(): just need to build and run (finite objf, cor
    # estimated off the boundary) -- a few outer iterations is enough
    .mT <- function() {
      ini({
        tka <- log(1); tcl <- log(0.3); tv <- log(7)
        eta.cl ~ 0.1; add.sd <- 0.6; ar1.cor <- 0.3; tdf <- fix(5)
      })
      model({
        ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp ~ add(add.sd) + ar(ar1.cor) + dt(tdf)
      })
    }
    .fT <- .nlmixr(.mT, .dat, est = "focei",
                   control = foceiControl(print = 0, maxOuterIterations = 10L))
    expect_true(is.finite(.fT$objf))
    expect_lt(unname(.fT$parFixedDf["ar1.cor", "Estimate"]), 0.95)

    .mC <- function() {
      ini({
        tka <- log(1); tcl <- log(0.3); tv <- log(7)
        eta.cl ~ 0.1; add.sd <- 0.6; ar1.cor <- 0.3
      })
      model({
        ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp ~ add(add.sd) + ar(ar1.cor) + dcauchy()
      })
    }
    .fC <- .nlmixr(.mC, .dat, est = "focei",
                   control = foceiControl(print = 0, maxOuterIterations = 10L))
    expect_true(is.finite(.fC$objf))
    expect_lt(unname(.fC$parFixedDf["ar1.cor", "Estimate"]), 0.95)
  })
})
