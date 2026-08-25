## The focei theta-reset path saves its working buffers into an environment and
## copies them back when it restarts (src/inner.cpp, saveIntoEnvrionment() /
## restoreFromEnvrionment()).  Each buffer's saved length comes from the length
## recorded at allocation, so a save covers exactly what was allocated: reading
## past the end of the eta block, or restoring a truncated buffer over live
## state, both corrupt the fit and can abort the process outright.  This file
## fits the scenario that resets its thetas; under an ASAN build it is what
## catches an out-of-bounds save directly.
nmTest({

  .thetaResetModel <- function() {
    ini({
      tka <- 0.45
      tcl <- 1.0
      tv <- 3.45
      eta.ka ~ 0.09
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  }

  # simulated from parameters well away from the initial estimates, which is
  # what walks the optimizer far enough to reset
  .thetaResetData <- function(nid = 6, seed = 1234) {
    rxode2::rxWithSeed(seed, {
      .ev <- rxode2::et(amt = 320, cmt = "depot", id = seq_len(nid))
      .ev <- rxode2::et(.ev, seq(0.5, 24, by = 1.5))
      .sim <- rxode2::rxSolve(.thetaResetModel(), .ev,
                              params = c(tka = 0.6, tcl = 1.1, tv = 3.6))
      .dat <- as.data.frame(.sim)[, c("id", "time", "cp")]
      .dat$cp <- .dat$cp + stats::rnorm(nrow(.dat), 0, 0.3)
      names(.dat) <- c("ID", "TIME", "DV")
      .dat$AMT <- 0
      .dat$EVID <- 0
      .dose <- data.frame(ID = seq_len(nid), TIME = 0, DV = NA, AMT = 320,
                          EVID = 1)
      .dat <- rbind(.dose, .dat)
      .dat[order(.dat$ID, .dat$TIME, -.dat$EVID), ]
    })
  }

  test_that("a focei fit that saves and restores its buffers completes", {
    .dat <- .thetaResetData()
    .ctl <- foceiControl(print = 0, covMethod = "analytic", sigdig = 6)
    .f1 <- .nlmixr(.thetaResetModel(), .dat, est = "focei", control = .ctl)
    .f2 <- .nlmixr(.thetaResetModel(), .dat, est = "focei", control = .ctl)

    expect_true(is.finite(.f1$objf))
    expect_true(all(is.finite(fixef(.f1))))
    expect_true(is.finite(.f2$objf))
    expect_true(all(is.finite(fixef(.f2))))
    # a restore that loses or overruns part of the saved state shows up as the
    # two fits disagreeing, so assert them equal rather than merely finite
    expect_equal(.f1$objf, .f2$objf)
    expect_equal(unname(fixef(.f1)), unname(fixef(.f2)))
  })

})
