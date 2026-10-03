## The focei theta-reset path saves its working buffers into an environment and
## copies them back when it restarts (src/inner.cpp, saveIntoEnvironment() /
## restoreFromEnvironment()).  Each buffer's length used to be re-derived from a
## copy of the allocation's formula, and every copy had drifted: three saved
## short, and the eta block's copy claimed more than had been allocated -- so
## saving read past the end of that block and restoring wrote the overshoot
## back, corrupting whatever followed it and, intermittently, aborting the R
## process.

nmTest({
  # A model whose CL is initialized three orders of magnitude off scale, so the
  # etas drift and the reset actually fires.  A fit that converges cleanly never
  # resets at all -- the reset is a recovery path -- so exercising it requires a
  # struggling fit.  Where this one lands is not the point and is not asserted;
  # what is asserted is that the reset ran, that the fit came back intact, and
  # that repeating it gives the same answer.
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
  # resetThetaP = 0 leaves resetThetaSize infinite, and thetaReset() returns on
  # std::isinf(size) before it can fire).  So it has to be switched on
  # explicitly AND counted -- without the count these tests would pass without
  # ever reaching the code they cover.
  #
  # NB: nlmixr() directly, not the .nlmixr() helper -- that helper wraps the fit
  # in suppressMessages(), whose handler muffles the reset messages before an
  # outer handler can count them.
  .fitCountingResets <- function(est, control) {
    .n <- 0L
    .fit <- withCallingHandlers(
      suppressWarnings(nlmixr(.driftingModel(), nlmixr2data::theo_sd, est = est, control = control)),
      message = function(m) {
        if (grepl("ETA drift", conditionMessage(m), fixed = TRUE)) {
          .n <<- .n + 1L
        }
        tryInvokeRestart("muffleMessage")
      }
    )
    list(fit = .fit, nReset = .n)
  }

  # Reaching the reset depends on the optimizer struggling.  It fires in every
  # configuration measured here, but a platform whose BLAS walks a different
  # trajectory could in principle converge straight through.  Skip loudly in
  # that case rather than fail -- and never pass quietly, which is what this
  # file did before.
  .skipIfNoReset <- function(...) {
    if (all(c(...) == 0L)) {
      skip("no theta reset was triggered here, so the save/restore path was never reached")
    }
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

  # Two identical runs must land in the same place.  This is what actually pins
  # the round trip: a save that drops the tail of a block, or a restore that
  # writes past the end of one, leaves state behind that differs between runs,
  # so the fits diverge grossly even though each on its own still looks finite.
  #
  # The tolerance is deliberate, and measured.  Repeating an identical
  # NON-resetting fit reproduces bit for bit, but a reset-heavy fit still drifts
  # by ~1e-6 relative between runs in the same session.  That drift is a
  # separate, older problem that this fix does not address, so it is tolerated
  # here rather than asserted away -- while still being four orders of magnitude
  # tighter than anything a truncated or overrunning restore would produce.
  .expectSameFit <- function(a, b) {
    expect_equal(a$objf, b$objf, tolerance = 1e-3)
    expect_equal(unname(fixef(a)), unname(fixef(b)), tolerance = 1e-3)
  }

  test_that("a focei fit that saves and restores its buffers completes", {
    .ctl <- foceiControl(
      resetThetaP = 0.2,
      resetThetaCheckPer = 1,
      print = 0,
      maxOuterIterations = 20L,
      covMethod = "",
      calcTables = FALSE
    )
    .a <- .fitCountingResets("focei", .ctl)
    .b <- .fitCountingResets("focei", .ctl)

    # the save/restore path ran
    .skipIfNoReset(.a$nReset, .b$nReset)
    expect_gt(.a$nReset, 0L)
    expect_gt(.b$nReset, 0L)
    # ... the fits came back intact, rather than aborting the process.  Before
    # the buffer lengths came from the allocation, this is where saving read
    # past the end of the eta block and restoring wrote the overshoot back;
    # that overwrite took down test-matexp.R outright, and an ASAN build flags
    # the read here directly.
    .expectWellFormed(.a$fit)
    .expectWellFormed(.b$fit)
    .expectSameFit(.a$fit, .b$fit)
  })

  test_that("a theta reset under nAGQ>0 restores its fit state", {
    # The fullTheta block is [fullTheta|theta|initPar|scaleC|aqx|aqw].  The
    # aqx/aqw tail is the nAGQ^neta Gauss-Hermite node grid, which setupAq1_()
    # refills from the fit environment on every entry -- before the restore runs
    # -- so it is deliberately left out of the save.  Pin that a reset under
    # nAGQ>0 still restores the part that matters, reproducibly.
    .ctl <- agqControl(
      nAGQ = 3,
      resetThetaP = 0.2,
      resetThetaCheckPer = 1,
      print = 0,
      maxOuterIterations = 20L,
      covMethod = "",
      calcTables = FALSE
    )
    .a <- .fitCountingResets("agq", .ctl)
    .b <- .fitCountingResets("agq", .ctl)

    .skipIfNoReset(.a$nReset, .b$nReset)
    expect_gt(.a$nReset, 0L)
    expect_gt(.b$nReset, 0L)
    .expectWellFormed(.a$fit)
    .expectWellFormed(.b$fit)
    .expectSameFit(.a$fit, .b$fit)
  })

  # The drifting model's first reset fires at the initial estimates.  With the
  # restart loop off it ends the fit, and .thetaReset keeps what it reset to.
  test_that("a theta reset moves each eta's drift into the theta it refers to", {
    skip_on_cran()
    # the etas at the initial estimates (the inner problem of a posthoc fit)
    .eta0 <- suppressMessages(suppressWarnings(nlmixr(
      .driftingModel(),
      nlmixr2data::theo_sd,
      est = "focei",
      control = foceiControl(print = 0, maxOuterIterations = 0L, covMethod = "", calcTables = FALSE)
    )))$eta
    .eta0 <- as.matrix(.eta0[, c("eta.ka", "eta.cl")])
    withr::local_options(nlmixr2.retryFocei = FALSE)
    rm(list = intersect(ls(.thetaReset), c("thetaIni", "omegaTheta", "etaMat", "nF")), envir = .thetaReset)
    expect_error(
      suppressMessages(suppressWarnings(nlmixr(
        .driftingModel(),
        nlmixr2data::theo_sd,
        est = "focei",
        control = foceiControl(
          resetThetaP = 0.2,
          resetThetaCheckPer = 1,
          print = 0,
          maxOuterIterations = 20L,
          covMethod = "",
          calcTables = FALSE
        )
      ))),
      "theta reset"
    )
    expect_equal(.thetaReset$nF, 0)
    # tv and add.sd have no eta: they stay at their initial estimates, as does
    # omega (chol(omega^-1) with sqrt diagonals: omega^(-1/4))
    expect_equal(.thetaReset$thetaIni[3:4], c(-1, 0.7))
    expect_equal(.thetaReset$omegaTheta, c(0.6, 0.3)^-0.25)
    # every subject keeps its ka and cl: each eta's drift is taken out of the
    # eta and put into tka and tcl, its own theta
    expect_equal(
      sweep(.thetaReset$etaMat, 2, .thetaReset$thetaIni[1:2], "+"),
      sweep(.eta0, 2, c(0.45, -3.2), "+"),
      tolerance = 1e-6,
      ignore_attr = TRUE
    )
    expect_true(all(abs(.thetaReset$thetaIni[1:2] - c(0.45, -3.2)) > 1))
  })

  test_that("the restart after a theta reset starts from the values it reset to", {
    skip_on_cran()
    rm(list = intersect(ls(.thetaReset), c("thetaIni", "omegaTheta", "etaMat", "nF")), envir = .thetaReset)
    # the final reset fires after the optimization, so its restart continues
    # the parameter history
    .f <- suppressMessages(suppressWarnings(nlmixr(
      .driftingModel(),
      nlmixr2data::theo_sd,
      est = "focei",
      control = foceiControl(
        resetThetaFinalP = 0.2,
        print = 0,
        maxOuterIterations = 20L,
        covMethod = "",
        calcTables = FALSE
      )
    )))
    if (is.null(.thetaReset$nF)) {
      skip("no theta reset was triggered here, so the restart was never reached")
    }
    expect_gt(.thetaReset$nF, 0)
    .un <- .f$parHistData[.f$parHistData$type == "Unscaled", ]
    .first <- unlist(.un[.un$iter == .thetaReset$nF + 1, c("tka", "tcl", "tv", "add.sd", "o1", "o2")])
    expect_equal(unname(.first), c(.thetaReset$thetaIni, .thetaReset$omegaTheta), tolerance = 1e-10)
    expect_true(is.finite(.f$objf))
  })
})
