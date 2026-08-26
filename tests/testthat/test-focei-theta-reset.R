## The focei theta-reset path saves its working buffers into an environment and
## copies them back when it restarts (src/inner.cpp, saveIntoEnvironment() /
## restoreFromEnvironment()).  Each buffer's length used to be re-derived from a
## copy of the allocation's formula, and every copy had drifted: three saved
## short, and the eta block's copy claimed more than had been allocated -- so
## saving read past the end of that block and restoring wrote the overshoot
## back, corrupting whatever followed it and, intermittently, aborting the R
## process.

nmTest({

  # A model whose CL is initialized well off scale, so the etas drift and the
  # reset actually fires.  A fit that converges cleanly never resets at all --
  # the reset is a recovery path, so exercising it requires a struggling fit.
  # The estimates it lands on are NOT the point here and are deliberately not
  # asserted; what is asserted is that the reset machinery ran and the fit came
  # back intact instead of corrupting the heap.
  .driftingModel <- function() {
    ini({
      tka <- 0.45
      tcl <- -3.2
      tv <- -1
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  # The buffers are only saved and restored when a fit actually resets its
  # thetas, and the ETA-drift reset is OFF by default (foceiControl's
  # resetThetaP = 0 leaves resetThetaSize infinite, and thetaReset() returns
  # before it can fire).  So it has to be switched on explicitly AND counted --
  # without the count these tests would pass without ever reaching the code
  # they cover.
  # NB: nlmixr() directly, not the .nlmixr() test helper -- that helper wraps the
  # fit in suppressMessages(), whose own handler muffles the reset messages
  # before an outer handler can count them.
  .fitCountingResets <- function(est, control) {
    .n <- 0L
    .fit <- withCallingHandlers(
      suppressWarnings(nlmixr(.driftingModel(), nlmixr2data::theo_sd,
                              est = est, control = control)),
      message = function(m) {
        if (grepl("ETA drift", conditionMessage(m), fixed = TRUE)) {
          .n <<- .n + 1L
        }
        tryInvokeRestart("muffleMessage")
      })
    list(fit = .fit, nReset = .n)
  }

  # A truncated or misaligned restore leaves fullTheta/theta/initPar/scaleC out
  # of step with each other, which shows up as a fit whose reported parameters
  # no longer match the model.
  .expectWellFormed <- function(fit) {
    expect_equal(names(fixef(fit)), c("tka", "tcl", "tv", "add.sd"))
    expect_true(all(is.finite(fixef(fit))))
    expect_true(is.finite(fit$objf))
    expect_true(all(is.finite(fit$omega)))
    expect_equal(dim(fit$omega), c(2L, 2L))
  }

  test_that("a focei fit that saves and restores its buffers completes", {
    .r <- .fitCountingResets(
      "focei",
      foceiControl(resetThetaP = 0.2, resetThetaCheckPer = 1, print = 0,
                   maxOuterIterations = 20L, covMethod = "",
                   calcTables = FALSE))

    # the save/restore path ran
    expect_gt(.r$nReset, 0L)
    # ... and the fit came back intact, rather than aborting the process.
    # Before the buffer lengths came from the allocation, this is where saving
    # read past the end of the eta block and restoring wrote the overshoot
    # back; that overwrite took down test-matexp.R outright, and an ASAN build
    # flags the read here directly.
    .expectWellFormed(.r$fit)
  })

  test_that("a theta reset under nAGQ>0 restores its fit state", {
    # The fullTheta block is [fullTheta|theta|initPar|scaleC|aqx|aqw].  The
    # aqx/aqw tail is the nAGQ^neta Gauss-Hermite node grid, which setupAq1_()
    # refills from the fit environment on every entry -- before the restore
    # runs -- so it is deliberately left out of the save.  Pin that a reset
    # under nAGQ>0 still restores the part that matters.
    .r <- .fitCountingResets(
      "agq",
      agqControl(nAGQ = 3, resetThetaP = 0.2, resetThetaCheckPer = 1,
                 print = 0, maxOuterIterations = 20L, covMethod = "",
                 calcTables = FALSE))

    expect_gt(.r$nReset, 0L)
    .expectWellFormed(.r$fit)
  })

})
