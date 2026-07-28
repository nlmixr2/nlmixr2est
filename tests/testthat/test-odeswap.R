# Shared solve-pool registry (src/odeSwap.cpp): the pool is sized for the
# largest-neq registered model, and a private lhs buffer is required exactly when
# the widest-lhs model is NOT that pool model.
#
# The pure .odeSwapPlanFor() checks need no fit, and one plain focei fit checks
# the registry against what rxode2 reports.  Everything fit-heavy lives in
# test-odeswap-fit.R (weekly): compiling several more models here evicts entries
# from rxode2's model cache and breaks code generation in a LATER test file --
# it showed up as "user function 'expit' failed to produce code" in test-nlm.R.

nmTest({
  test_that("the pool plan picks max neq, with a deterministic tie-break", {
    # single model: it is the pool, nothing to override, no scratch
    p <- .odeSwapPlanFor(c(4L), c(6L))
    expect_identical(p$poolSlot, 0L)
    expect_identical(p$poolNeq, 4L)
    expect_identical(p$nLoaded, 1L)
    expect_false(p$overrideNeeded)
    expect_identical(p$scratchNlhs, 0L)

    # largest neq wins regardless of position, and a smaller peer needs an override
    p <- .odeSwapPlanFor(c(4L, 26L, 2L), c(6L, 29L, 4L))
    expect_identical(p$poolSlot, 1L)
    expect_identical(p$poolNeq, 26L)
    expect_identical(p$poolNlhs, 29L)
    expect_true(p$overrideNeeded)
    expect_identical(p$nLoaded, 3L)

    # neq tie -> the wider lhs wins, which is what keeps scratchNlhs at 0
    p <- .odeSwapPlanFor(c(10L, 10L), c(4L, 20L))
    expect_identical(p$poolSlot, 1L)
    expect_identical(p$scratchNlhs, 0L)
    expect_false(p$overrideNeeded)

    # fully tied -> lowest slot, so the choice never depends on source order
    p <- .odeSwapPlanFor(c(10L, 10L), c(7L, 7L))
    expect_identical(p$poolSlot, 0L)

    # an all-zero slot is unloaded and ignored; a 0-state slot WITH lhs is not
    p <- .odeSwapPlanFor(c(0L, 0L, 5L), c(0L, 99L, 3L))
    expect_identical(p$poolSlot, 2L)     # only slot 2 has states
    expect_identical(p$nLoaded, 2L)      # slots 1 and 2
    expect_identical(p$maxNlhs, 99L)     # slot 1 is the widest
    expect_identical(p$scratchNlhs, 99L)
    expect_true(p$overrideNeeded)        # slot 1 must compact from 5 to 0

    # nothing loaded at all
    p <- .odeSwapPlanFor(c(0L, 0L), c(0L, 0L))
    expect_identical(p$poolSlot, -1L)
    expect_identical(p$nLoaded, 0L)
    expect_identical(p$scratchNlhs, 0L)
  })

  test_that("a loaded model with zero ODE states still counts", {
    # A solved-form (linCmt) or purely algebraic model has neq == 0 but real lhs
    # outputs.  Treating neq == 0 as "not loaded" would drop it from maxNlhs and
    # silently skip the scratch buffer its calc_lhs needs.
    p <- .odeSwapPlanFor(c(2L, 0L), c(4L, 15L))
    expect_identical(p$nLoaded, 2L)
    expect_identical(p$poolSlot, 0L)
    expect_identical(p$maxNlhs, 15L)
    expect_identical(p$scratchNlhs, 15L)
    expect_true(p$needsScratch)
    expect_true(p$overrideNeeded)   # the 0-state model must compact the stride

    # and a zero-state model is still a valid pool when it is all there is
    p <- .odeSwapPlanFor(c(0L), c(9L))
    expect_identical(p$poolSlot, 0L)
    expect_identical(p$poolNeq, 0L)
    expect_identical(p$scratchNlhs, 0L)
  })

  test_that("a wider-lhs peer that is not the pool model forces the scratch buffer", {
    # The load-bearing case: rxode2's per-thread ind->lhs slice is exactly
    # op->nlhs (= poolNlhs) wide and rxode2 has no lhs override, so reading the
    # wider model through that slice would run past this thread's region.
    # est="impmap" reaches this with an ordinary combined-error model: residual
    # thetas add d(V)/d(theta) lhs columns but no sensitivity states.
    p <- .odeSwapPlanFor(c(10L, 3L), c(4L, 20L))
    expect_identical(p$poolSlot, 0L)      # 10 states sizes the pool
    expect_identical(p$maxNlhsSlot, 1L)   # but slot 1 has the widest lhs
    expect_identical(p$maxNlhs, 20L)
    expect_identical(p$scratchNlhs, 20L)  # so the read needs our own buffer
    expect_true(p$overrideNeeded)

    # and it is not needed when the pool model is also the widest
    p <- .odeSwapPlanFor(c(10L, 3L), c(20L, 4L))
    expect_identical(p$scratchNlhs, 0L)
  })

  test_that("three augmented models coexist: pool by max ODEs, then the lhs pointer", {
    # The analytic path compiles up to three augmented models for one fit (order-2
    # gradient, order-1 AGQ node, covariance over its own direction set).  They all
    # stay registered together; solving one calls ITS entry points.  The pool is
    # built once for whichever has the most ODE states, and only then is the lhs
    # pointer decided.
    # slots: inner, pred, <unused>, <unused>, outer, outerNode, outerCov
    #   outer     26 states / 29 lhs   <- most ODEs, so it sizes the pool
    #   outerNode 14 / 17
    #   outerCov  20 / 23
    p <- .odeSwapPlanFor(c(8L, 2L, 0L, 0L, 26L, 14L, 20L),
                         c(8L, 2L, 0L, 0L, 29L, 17L, 23L))
    expect_identical(p$nLoaded, 5L)
    expect_identical(p$poolSlot, 4L)      # the order-2 gradient model
    expect_identical(p$poolNeq, 26L)
    expect_true(p$overrideNeeded)         # every other model runs compacted
    # widest lhs belongs to the pool model, so rxode2's own slice is wide enough
    expect_identical(p$maxNlhsSlot, 4L)
    expect_identical(p$scratchNlhs, 0L)
    expect_false(p$needsScratch)

    # same three models, but the covariance model carries extra outputs without
    # extra states -- now the widest lhs is NOT the pool model and reads of it
    # need a private buffer
    p <- .odeSwapPlanFor(c(8L, 2L, 0L, 0L, 26L, 14L, 20L),
                         c(8L, 2L, 0L, 0L, 29L, 17L, 40L))
    expect_identical(p$poolSlot, 4L)      # pool choice is unchanged: still max ODEs
    expect_identical(p$maxNlhsSlot, 6L)   # but the cov model is the widest reader
    expect_identical(p$scratchNlhs, 40L)
    expect_true(p$needsScratch)
  })

  test_that("the registry reports each peer's true neq/nlhs for a plain focei fit", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    f <- suppressMessages(nlmixr2(one.cmt, nlmixr2data::theo_sd, "focei",
                                  foceiControl(print = 0L, covMethod = "",
                                               maxOuterIterations = 0L,
                                               calcTables = FALSE)))
    i <- .odeSwapInfo()
    m <- i$models
    expect_true(is.data.frame(m))
    expect_identical(nrow(m), length(.odeSwapSlots))

    # exactly one model sizes the pool, and it is the largest registered one
    expect_identical(sum(m$sizesPool), 1L)
    expect_identical(i$poolNeq, max(m$neq))

    # the registry's counts must match what rxode2 reports for the same models,
    # read independently from R.  Unloaded slots carry NA names, so match with
    # %in% (== would propagate NA into the subscript).
    mods <- f$foceiModel
    expect_identical(m$neq[m$name %in% "inner"],
                     length(rxode2::rxModelVars(mods$inner)$state))
    expect_identical(m$nlhs[m$name %in% "inner"],
                     length(rxode2::rxModelVars(mods$inner)$lhs))
    expect_identical(m$neq[m$name %in% "pred"],
                     length(rxode2::rxModelVars(mods$predNoLhs)$state))

    # first use of getOpNlhs() in this package -- assert it, do not assume it
    expect_identical(i$opNeq, i$poolNeq)
    expect_identical(i$opNlhs, i$poolNlhs)

    # a plain focei fit: the inner model is the pool and nothing is wider
    expect_identical(i$poolName, "inner")
    expect_false(i$needsScratch)
  })
})
