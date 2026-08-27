## Fitting the same data twice in one session has to give the same answer
## (nlmixr2est issue 1020).
##
## FOCEi solves one subject at a time through rxode2's `ind_solve()` and
## re-sorts the subjects around every inner-optimization and every gradient
## loop.  Two rxode2 bugs met there:
##
##   * `sortIds()` ordered the subjects most-expensive-first from
##     `ind->solveTime`, the accumulated `clock()` time of their earlier
##     solves, so the order was a function of wall-clock timing -- it differed
##     between two runs of the same problem and between successive passes of
##     one run.
##   * `ind_solve()`'s argument is a subject id, but several of rxode2's
##     per-individual drivers read it as a position in `rx->ordId`.  While
##     `ordId` is the identity the two agree; once it is not, the wrong
##     individual is integrated.
##
## The reorder only fires when there are at least `throttle` (2) times more
## threads than subjects, so past that threshold the likelihood, its gradient
## and hence the whole fit followed a solve order that came from thread
## timing.  Both are fixed in rxode2; these tests are the gate.
##
## Two things the tests have to do deliberately, or they prove nothing:
##
##   * raise the thread count to 2x the subject count for their own duration,
##     because the reorder does not trigger below it, and
##   * give the subjects DIFFERENT numbers of samples.  Equal-sized subjects
##     sort to the identity order, and an identity `ordId` hides the
##     id-vs-position confusion completely.

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

  .ini <- c(tka = 0.45, tcl = 1.0, tv = 3.45, add.sd = 0.7)

  # simulated from parameters away from the initial estimates, so the optimizer
  # has to travel -- a fit that stops where it started hides the problem.
  # Subject i gets 8 + 2i samples: the sizes have to differ for the solve order
  # to be anything but the identity.
  .reproData <- function(model, nid, seed = 1234) {
    rxode2::rxWithSeed(seed, {
      .ev <- do.call(
        rbind,
        lapply(seq_len(nid), function(.i) {
          rbind(data.frame(id = .i, time = 0, amt = 320, evid = 1),
                data.frame(id = .i,
                           time = seq(0.5, 24, length.out = 8L + 2L * .i),
                           amt = 0, evid = 0))
        }))
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

  # How many subjects this machine can exercise the reorder with: it needs
  # `throttle` (2) threads per subject, so the widest useful problem is half
  # the logical CPUs.  Capped at `want`, which is a big enough problem for the
  # estimates to mean something.
  .reproNid <- function(want = 6L) {
    .old <- rxode2::getRxThreads()
    on.exit(rxode2::setRxThreads(.old), add = TRUE)
    rxode2::setRxThreads(0L) # 0 means every logical CPU
    max(2L, min(want, rxode2::getRxThreads() %/% 2L))
  }

  # Fit twice at the raised thread count and once at two threads, and demand
  # all three agree.  `skip_if()` is on what `setRxThreads()` actually granted
  # -- it clamps to the machine silently, and `getRxThreads()` reports the
  # current SETTING rather than the number of CPUs, so asking before setting
  # skips machines that could have run this.
  .expectReproducible <- function(model, nid) {
    .threads <- 2L * nid
    .old <- rxode2::getRxThreads()
    on.exit(rxode2::setRxThreads(.old), add = TRUE)
    rxode2::setRxThreads(.threads)
    skip_if(rxode2::getRxThreads() < .threads,
            paste0("needs ", .threads, " threads (2 per subject) to reorder"))

    .dat <- .reproData(model, nid)
    .ctl <- foceiControl(print = 0, covMethod = "")
    .f1 <- .nlmixr(model, .dat, est = "focei", control = .ctl)
    .f2 <- .nlmixr(model, .dat, est = "focei", control = .ctl)

    expect_equal(.f1$objf, .f2$objf)
    expect_equal(unname(fixef(.f1)), unname(fixef(.f2)))

    # ... and the answer must not depend on how many threads are available
    rxode2::setRxThreads(2L)
    .f3 <- .nlmixr(model, .dat, est = "focei", control = .ctl)
    expect_equal(.f1$objf, .f3$objf)
    expect_equal(unname(fixef(.f1)), unname(fixef(.f3)))

    # the fit also has to have gone somewhere: the failure mode was an
    # optimizer that gave up at (or beside) the initial estimates
    expect_false(isTRUE(all.equal(unname(fixef(.f1)[names(.ini)]),
                                  unname(.ini), tolerance = 1e-4)))
    # ... and somewhere right.  Only worth asking of the full-size problem;
    # a narrow machine runs this with too few subjects to pin the estimates.
    if (nid >= 6L) {
      expect_equal(unname(fixef(.f1)[c("tka", "tcl", "tv")]),
                   c(0.6, 1.1, 3.6), tolerance = 0.1)
    }
    invisible(.f1)
  }

  test_that("the same focei fit twice in one session gives the same answer", {
    .expectReproducible(.reproModel(), .reproNid())
  })

  test_that("a linCmt() focei fit twice in one session gives the same answer", {
    .expectReproducible(.reproModelLinCmt(), .reproNid())
  })

})
