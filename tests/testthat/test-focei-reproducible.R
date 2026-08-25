## Fitting the same data twice in one session has to give the same answer
## (nlmixr2est issue 1020).
##
## FOCEi solves one subject at a time through rxode2's `ind_solve()`, whose
## argument is a subject id.  Several of rxode2's per-individual drivers read
## that argument as a position in `rx->ordId` instead, and `rx->ordId` stops
## being the identity once rxode2 reorders subjects most-expensive-first --
## which it does when there are at least `throttle` (2) times more threads than
## subjects.  Past that threshold the wrong individual was integrated, so the
## likelihood, its gradient, and hence the whole fit depended on a solve order
## that comes from wall-clock timing.  Fixed in rxode2; these tests need that
## fix to pass.
##
## The tests need threads >= 2 * subjects to exercise the path at all, so they
## raise the thread count for their own duration.

nmTest({

  .reproModel <- function() {
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

  # the same structural model as an analytic solution, which reaches a
  # different family of rxode2 drivers (ind_linCmt0 / ind_linCmt0H)
  .reproModelLinCmt <- function() {
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
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }

  # simulated from parameters away from the initial estimates, so the optimizer
  # has to travel -- a fit that stops where it started hides the problem
  .reproData <- function(model, nid = 6, seed = 1234) {
    rxode2::rxWithSeed(seed, {
      .ev <- rxode2::et(amt = 320, cmt = "depot", id = seq_len(nid))
      .ev <- rxode2::et(.ev, seq(0.5, 24, by = 1.5))
      .sim <- rxode2::rxSolve(model, .ev,
                              params = c(tka = 0.6, tcl = 1.1, tv = 3.6))
      .dat <- as.data.frame(.sim)[, c("id", "time", "cp")]
      .dat$cp <- .dat$cp + stats::rnorm(nrow(.dat), 0, 0.3)
      names(.dat) <- c("ID", "TIME", "DV")
      .dat$AMT <- 0
      .dat$EVID <- 0
      .dat <- rbind(data.frame(ID = seq_len(nid), TIME = 0, DV = NA,
                               AMT = 320, EVID = 1),
                    .dat)
      .dat[order(.dat$ID, .dat$TIME, -.dat$EVID), ]
    })
  }

  .ini <- c(tka = 0.45, tcl = 1.0, tv = 3.45, add.sd = 0.7)

  test_that("the same focei fit twice in one session gives the same answer", {
    .nid <- 6L
    .dat <- .reproData(.reproModel(), .nid)
    # the reorder only triggers at cores >= subjects * throttle (throttle is 2)
    .threads <- 2L * .nid
    skip_if(rxode2::rxCores() < .threads,
            "needs at least 2 threads per subject to exercise the reorder")
    .old <- rxode2::getRxThreads()
    on.exit(rxode2::setRxThreads(.old))
    rxode2::setRxThreads(.threads)

    .ctl <- foceiControl(print = 0, covMethod = "")
    .f1 <- .nlmixr(.reproModel(), .dat, est = "focei", control = .ctl)
    .f2 <- .nlmixr(.reproModel(), .dat, est = "focei", control = .ctl)

    expect_equal(.f1$objf, .f2$objf)
    expect_equal(unname(fixef(.f1)), unname(fixef(.f2)))

    # ... and the answer must not depend on how many threads are available
    rxode2::setRxThreads(2L)
    .f3 <- .nlmixr(.reproModel(), .dat, est = "focei", control = .ctl)
    expect_equal(.f1$objf, .f3$objf)
    expect_equal(unname(fixef(.f1)), unname(fixef(.f3)))

    # the fit also has to have gone somewhere: the failure mode was an
    # optimizer that gave up at (or beside) the initial estimates
    expect_false(isTRUE(all.equal(unname(fixef(.f1)[names(.ini)]),
                                  unname(.ini), tolerance = 1e-4)))
    # it converges to what the data was simulated from
    expect_equal(unname(fixef(.f1)[c("tka", "tcl", "tv")]),
                 c(0.6, 1.1, 3.6), tolerance = 0.1)
  })

  test_that("a linCmt() focei fit twice in one session gives the same answer", {
    .nid <- 6L
    .dat <- .reproData(.reproModelLinCmt(), .nid)
    .threads <- 2L * .nid
    skip_if(rxode2::rxCores() < .threads,
            "needs at least 2 threads per subject to exercise the reorder")
    .old <- rxode2::getRxThreads()
    on.exit(rxode2::setRxThreads(.old))
    rxode2::setRxThreads(.threads)

    .ctl <- foceiControl(print = 0, covMethod = "")
    .f1 <- .nlmixr(.reproModelLinCmt(), .dat, est = "focei", control = .ctl)
    .f2 <- .nlmixr(.reproModelLinCmt(), .dat, est = "focei", control = .ctl)

    expect_equal(.f1$objf, .f2$objf)
    expect_equal(unname(fixef(.f1)), unname(fixef(.f2)))

    rxode2::setRxThreads(2L)
    .f3 <- .nlmixr(.reproModelLinCmt(), .dat, est = "focei", control = .ctl)
    expect_equal(.f1$objf, .f3$objf)
    expect_equal(unname(fixef(.f1)), unname(fixef(.f3)))

    expect_false(isTRUE(all.equal(unname(fixef(.f1)[names(.ini)]),
                                  unname(.ini), tolerance = 1e-4)))
    expect_equal(unname(fixef(.f1)[c("tka", "tcl", "tv")]),
                 c(0.6, 1.1, 3.6), tolerance = 0.1)
  })

})
