# Combined eta+theta sensitivity build (#958) applied to imp/impmap's M-step.
#
# est="imp"/"impmap" update their non-mu structural + residual-error thetas by
# a Newton step whose score/Hessian needs d(f)/d(theta) and d(V)/d(theta) at
# every E-step sample.  By default that is a SECOND, dedicated rxThetaSens
# model re-solving every sample the E-step had just solved through the inner
# model to get the importance weights -- two ODE integrations per sample.
# impmapControl(combSens=TRUE) (opt-in, default FALSE) instead carries the
# theta columns on the inner model itself, and impEStep (src/imp.cpp)
# harvests them straight off its own per-sample inner solve
# (impThetaSensCollect's reuseSolve path), so the common case (sir=FALSE)
# needs only ONE solve per sample.
#
# combSens=TRUE is opt-in, not the default: folding the theta-sensitivity
# states into the SAME coupled ODE system as the eta sensitivities changes the
# integrator's adaptive step-size path, so results move by a small amount
# relative to the two-model default -- not a correctness issue (see the
# FOCEI-convergence check below), but enough to occasionally cross a
# diagnostic threshold, so it must not change any existing fit's numbers
# silently.
#
# These tests assert the MECHANISM is used (the harvest counter moves), not
# just that results look plausible -- a wrong wiring could still converge to
# something.
nmTest({

  .harvestN <- function() nlmixr2est:::.odeSwapInfo()$impThetaSensHarvestN

  # tv is a non-mu STRUCTURAL theta (an ODE-state sensitivity, not just an
  # algebraic one); add.sd is a non-mu SIGMA theta.  Both get d(f)/d(theta) /
  # d(V)/d(theta) columns, exercising the full harvested M-step gradient.
  .mod <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("combSens defaults off: impThetaSensIdx alone does not request the fused build", {
    .d <- nlmixr2data::theo_sd
    expect_false(impmapControl()$combSens)
    .n0 <- .harvestN()
    rxode2::rxSetSeed(42)
    .f <- suppressWarnings(
      nlmixr2(.mod, .d, "imp", impControl(print = 0L, nIter = 5L, isample = 50L)))
    expect_true(inherits(.f, "nlmixr2FitCore"))
    # no fused build requested -> nothing harvested (unchanged default behavior)
    expect_equal(.harvestN(), .n0)
  })

  test_that("imp(combSens=TRUE) harvests theta-sensitivities from the E-step's own inner solve", {
    .d <- nlmixr2data::theo_sd
    .ff <- suppressWarnings(nlmixr2(.mod, .d, "focei", foceiControl(print = 0L, covMethod = "")))
    .n0 <- .harvestN()
    rxode2::rxSetSeed(42)
    .fh <- suppressWarnings(
      nlmixr2(.mod, .d, "imp",
              impControl(print = 0L, nIter = 30L, isample = 300L, combSens = TRUE)))
    expect_true(inherits(.fh, "nlmixr2FitCore"))
    # the mechanism actually ran (not just "results look fine")
    expect_true(.harvestN() > .n0)
    # combSens builds ONE fused model -- no separate thetaSens peer registered
    .info <- nlmixr2est:::.odeSwapInfo()
    expect_false("thetaSens" %in% .info$models$name)
    # the structural theta (tv) and sigma theta (add.sd) went through the
    # harvested M-step gradient and converged near FOCEI
    expect_equal(unname(fixef(.fh)["tv"]), unname(fixef(.ff)["tv"]), tolerance = 0.03)
    expect_equal(unname(fixef(.fh)["add.sd"]), unname(fixef(.ff)["add.sd"]), tolerance = 0.05)

    # sir=TRUE disables harvesting even with combSens=TRUE (impOuter's
    # harvestSens gate is off under any sir=TRUE), so the M-step falls back to
    # the un-fused solve-based read against the SAME fused inner model.  The
    # two paths read the same per-sample f/V/d(f)/d(theta)/d(V)/d(theta) off
    # the same ODE system, just at different points in the loop, so they
    # should land in the same neighborhood (not bit-identical: an EM run this
    # long accumulates its own floating-point path noise regardless of
    # harvesting, plus sir=TRUE's own RNG draws are not exactly the same
    # sequence).
    .n1 <- .harvestN()
    rxode2::rxSetSeed(42)
    .fs <- suppressWarnings(
      nlmixr2(.mod, .d, "imp",
              impControl(print = 0L, nIter = 30L, isample = 300L, combSens = TRUE,
                        sir = TRUE, sirSample = 300L)))
    expect_true(inherits(.fs, "nlmixr2FitCore"))
    # the fallback path harvested nothing new
    expect_equal(.harvestN(), .n1)
    expect_equal(unname(fixef(.fh)), unname(fixef(.fs)), tolerance = 0.02)
  })

  test_that("impmap(combSens=TRUE) (MAP search) harvests theta-sensitivities the same way", {
    .d <- nlmixr2data::theo_sd
    .n0 <- .harvestN()
    rxode2::rxSetSeed(42)
    .fh <- suppressWarnings(
      nlmixr2(.mod, .d, "impmap",
              impmapControl(print = 0L, nIter = 10L, isample = 50L, combSens = TRUE)))
    expect_true(inherits(.fh, "nlmixr2FitCore"))
    expect_true(.harvestN() > .n0)
  })

  test_that("no non-mu thetas: combSens(TRUE) has nothing to fuse, harvest counter unmoved", {
    # every theta here is mu-referenced (or fixed), so impmap has nothing to
    # differentiate -- combSens should resolve to off (nothing to fuse) even
    # though it was requested, and the harvest counter must not move.
    mmu <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- fix(3.45)
        eta.ka ~ 0.6; eta.cl ~ 0.3
        add.sd <- fix(0.7)
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .d <- nlmixr2data::theo_sd
    .n0 <- .harvestN()
    rxode2::rxSetSeed(42)
    .f <- suppressWarnings(
      nlmixr2(mmu, .d, "imp",
              impControl(print = 0L, nIter = 5L, isample = 50L, combSens = TRUE)))
    expect_true(inherits(.f, "nlmixr2FitCore"))
    expect_equal(.harvestN(), .n0)
  })

})
