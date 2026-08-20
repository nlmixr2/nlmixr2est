# Combined eta+theta sensitivity build (#958) applied to imp/impmap's M-step.
#
# est="imp"/"impmap" update their non-mu structural + residual-error thetas by
# a Newton step whose score/Hessian needs d(f)/d(theta) and d(V)/d(theta) at
# every E-step sample.  Without combSens that is a SECOND, dedicated
# rxThetaSens model re-solving every sample the E-step had just solved through
# the inner model to get the importance weights -- two ODE integrations per
# sample.  impmapControl(combSens=TRUE) (the default) instead carries the
# theta columns on the inner model itself, and impEStep (src/imp.cpp)
# harvests them straight off its own per-sample inner solve
# (impThetaSensCollect's reuseSolve path), so the common case (sir=FALSE)
# needs only ONE solve per sample.
#
# combSens=TRUE is the default because it agrees exactly (to floating-point
# noise) with the old two-model path -- see the ndiff-fix tests below for why
# that is worth pinning explicitly, not assumed: a linCmt() model with a
# SINGLE non-mu structural theta used to disagree measurably between the two
# builds, and that turned out to be an unrelated, pre-existing bug in
# src/odeSwap.cpp's peer-model swapping (rx->ndiff never restored per slot),
# not something combSens introduced.  impmapControl(combSens=FALSE) still
# works, for anyone who wants the two-model path.
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

  test_that("combSens defaults on: impThetaSensIdx alone requests the fused build", {
    .d <- nlmixr2data::theo_sd
    expect_true(impmapControl()$combSens)
    .n0 <- .harvestN()
    rxode2::rxSetSeed(42)
    .f <- suppressWarnings(
      nlmixr2(.mod, .d, "imp", impControl(print = 0L, nIter = 5L, isample = 50L)))
    expect_true(inherits(.f, "nlmixr2FitCore"))
    # fused build requested by default -> the harvest mechanism ran
    expect_true(.harvestN() > .n0)
    .info <- nlmixr2est:::.odeSwapInfo()
    expect_false("thetaSens" %in% .info$models$name)
  })

  test_that("imp harvests theta-sensitivities from the E-step's own inner solve", {
    .d <- nlmixr2data::theo_sd
    .ff <- suppressWarnings(nlmixr2(.mod, .d, "focei", foceiControl(print = 0L, covMethod = "")))
    .n0 <- .harvestN()
    rxode2::rxSetSeed(42)
    .fh <- suppressWarnings(
      nlmixr2(.mod, .d, "imp",
              impControl(print = 0L, nIter = 30L, isample = 300L)))
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

    # sir=TRUE disables harvesting (impOuter's harvestSens gate is off under
    # any sir=TRUE), so the M-step falls back to the un-fused solve-based read
    # against the SAME fused inner model.  The two paths read the same
    # per-sample f/V/d(f)/d(theta)/d(V)/d(theta) off the same ODE system, just
    # at different points in the loop, so they should land in the same
    # neighborhood (not bit-identical: an EM run this long accumulates its own
    # floating-point path noise regardless of harvesting, plus sir=TRUE's own
    # RNG draws are not exactly the same sequence).
    .n1 <- .harvestN()
    rxode2::rxSetSeed(42)
    .fs <- suppressWarnings(
      nlmixr2(.mod, .d, "imp",
              impControl(print = 0L, nIter = 30L, isample = 300L,
                        sir = TRUE, sirSample = 300L)))
    expect_true(inherits(.fs, "nlmixr2FitCore"))
    # the fallback path harvested nothing new
    expect_equal(.harvestN(), .n1)
    expect_equal(unname(fixef(.fh)), unname(fixef(.fs)), tolerance = 0.02)
  })

  test_that("impmap (MAP search) harvests theta-sensitivities the same way", {
    .d <- nlmixr2data::theo_sd
    .n0 <- .harvestN()
    rxode2::rxSetSeed(42)
    .fh <- suppressWarnings(
      nlmixr2(.mod, .d, "impmap",
              impmapControl(print = 0L, nIter = 10L, isample = 50L)))
    expect_true(inherits(.fh, "nlmixr2FitCore"))
    expect_true(.harvestN() > .n0)
  })

  test_that("combSens=FALSE still works and requests no fused build", {
    .d <- nlmixr2data::theo_sd
    .n0 <- .harvestN()
    rxode2::rxSetSeed(42)
    .f <- suppressWarnings(
      nlmixr2(.mod, .d, "imp",
              impControl(print = 0L, nIter = 5L, isample = 50L, combSens = FALSE)))
    expect_true(inherits(.f, "nlmixr2FitCore"))
    # opted out of the fused build -> nothing harvested
    expect_equal(.harvestN(), .n0)
  })

  test_that("no non-mu thetas: combSens has nothing to fuse, harvest counter unmoved", {
    # every theta here is mu-referenced (or fixed), so impmap has nothing to
    # differentiate -- combSens should resolve to off (nothing to fuse) even
    # though it defaults on, and the harvest counter must not move.
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
              impControl(print = 0L, nIter = 5L, isample = 50L)))
    expect_true(inherits(.f, "nlmixr2FitCore"))
    expect_equal(.harvestN(), .n0)
  })

  test_that("mixture (Nmix>1): harvested and solve-based M-steps agree on the expanded pseudo-subject mapping", {
    # impEStep/impOuter index a mixture's expanded pseudo-subjects as
    # id = i + j*nsub (j = 0-based component).  outSens is sized/filled by
    # that SAME id, so a wrong offset here would misalign a component's
    # samples against another component's zk weights -- this model has THREE
    # non-mu thetas that are not eta-linked (tka, p1, add.sd all lack a `+eta`
    # term), so d(f)/d(theta) actually differs across the two mixture
    # components' branches, not just a shared constant.
    .mkg <- function(cl0, ids) {
      ka <- 1.5; v <- 8
      do.call(rbind, lapply(ids, function(id) {
        cli <- cl0 * exp(stats::rnorm(1, 0, 0.2))
        tt <- c(0.5, 2, 6, 12)
        cp <- (100 * ka / (v * (ka - cli / v))) * (exp(-cli / v * tt) - exp(-ka * tt))
        cp <- pmax(cp, 1e-3) * exp(stats::rnorm(length(tt), 0, 0.1))
        rbind(data.frame(id = id, time = 0, dv = NA_real_, amt = 100, evid = 1, cmt = "depot"),
              data.frame(id = id, time = tt, dv = cp, amt = 0, evid = 0, cmt = "cen"))
      }))
    }
    .testSeed(11); rxode2::rxSetSeed(11)
    .d <- rbind(.mkg(3.0, 1:6), .mkg(9.0, 7:12))
    .d <- .d[order(.d$id, .d$time, -.d$evid), ]
    mmix <- function() {
      ini({
        tka <- log(1.5); tcl1 <- log(2.5); tcl2 <- log(8); tv <- log(8)
        p1 <- 0.5
        eta.cl ~ 0.2
        add.sd <- 0.3
      })
      model({
        ka <- exp(tka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(cen) <- ka * depot - cl / v * cen
        cp <- cen / v
        cp ~ add(add.sd)
      })
    }
    .n0 <- .harvestN()
    rxode2::rxSetSeed(42)
    .fh <- suppressWarnings(
      nlmixr2(mmix, .d, "impmap",
              impmapControl(print = 0L, nIter = 8L, isample = 40L)))
    expect_true(inherits(.fh, "nlmixr2FitCore"))
    expect_true(all(is.finite(fixef(.fh))))
    expect_true(.harvestN() > .n0)

    rxode2::rxSetSeed(42)
    .fs <- suppressWarnings(
      nlmixr2(mmix, .d, "impmap",
              impmapControl(print = 0L, nIter = 8L, isample = 40L, combSens = FALSE)))
    expect_true(inherits(.fs, "nlmixr2FitCore"))
    # A misaligned i+j*nsub mapping would scramble which component's samples
    # feed which component's Newton step -- the two components' clearances
    # (tcl1 low, tcl2 high) would no longer track between the harvested and
    # solve-based runs even though both start from the same seed/data.  Not
    # bit-identical (the two builds' ODE integrator paths differ slightly);
    # a loose tolerance is enough to catch a scrambled mapping, which would
    # miss by an order of magnitude or swap the two components.
    expect_equal(unname(fixef(.fh)[c("tcl1", "tcl2", "p1")]),
                unname(fixef(.fs)[c("tcl1", "tcl2", "p1")]), tolerance = 0.15)
  })

  test_that("odeSwapSolveInd restores rx->ndiff per peer (linCmtB Jacobian-cache fix)", {
    # Regression test for the bug combSens=TRUE originally uncovered on this
    # exact model shape: linCmt(), ONE non-mu structural theta (tka has no
    # eta), the OTHER two structural parameters (tv, ka's amplitude) pure
    # constants.  linCmtB()'s cached Jacobian (J/Jg) is keyed off rx->ndiff, a
    # single field on the shared rx_solve struct that used to be set ONCE
    # (whichever model first called rxSolve_()) and never restored when
    # odeSwap swapped to a DIFFERENT registered peer -- so impmap's inner
    # Hessian (impGetHessian -> calcEtaHessian) could read a stale,
    # iteration-old d(pred)/d(eta) left by a sibling peer with a different (or
    # absent) sensitivity need.  The practical symptom: impmap's proposal
    # covariance came out far too wide (a near-singular Hessian), which
    # LOOKED like an unusually healthy Pareto k-hat (over-coverage hides tail
    # problems) while actually wasting samples and masking a real tail issue.
    #
    # combSens=TRUE and combSens=FALSE are asserted to agree here -- before
    # the fix they did not (one read k-hat ~ -1.5, healthy; the other ~ 3,
    # unreliable, for the SAME subject).
    .pk <- function() {
      ini({tka <- 0.45; tcl <- 1; tv <- 3.45; eta.cl ~ 0.1; add.sd <- 0.7})
      model({ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
             linCmt() ~ add(add.sd)})
    }
    .d <- nlmixr2data::theo_sd
    run1 <- function(combSens) {
      rxode2::rxSetSeed(1)
      f <- suppressWarnings(nlmixr2(.pk, .d, "impmap",
                    impmapControl(print = 0L, nIter = 12L, isample = 300L,
                                  covMethod = "", auto = FALSE, gammaRule = "floor",
                                  combSens = combSens)))
      k <- f$env$impPsisK
      list(objf = as.numeric(f$objf), maxK = max(k, na.rm = TRUE),
           nAbove = sum(k > 0.7, na.rm = TRUE))
    }
    rF <- run1(FALSE)
    rT <- run1(TRUE)
    expect_equal(rF$objf, rT$objf, tolerance = 1e-4)
    expect_equal(rF$maxK, rT$maxK, tolerance = 1e-3)
    expect_equal(rF$nAbove, rT$nAbove)
    # the registry actually recorded a real (nonzero) ndiff for the inner
    # slot, AND a DIFFERENT (zero) one for the pred peer -- proof the fix has
    # something to restore between, not a vacuously-passing comparison on a
    # model where every registered slot needed the same ndiff anyway.  pred
    # is the doFD fallback slot (fInd->doFD, likInner0) -- a subject that
    # falls back to it mid-fit is exactly the scenario that used to leave
    # rx->ndiff on the wrong peer's value for whichever slot solved next.
    .info <- nlmixr2est:::.odeSwapInfo()
    .inner <- .info$models[.info$models$name %in% "inner", ]
    .pred <- .info$models[.info$models$name %in% "pred", ]
    expect_true(nrow(.inner) == 1L && isTRUE(.inner$ndiffSet) && .inner$ndiff > 0L)
    expect_true(nrow(.pred) == 1L && isTRUE(.pred$ndiffSet))
    expect_false(isTRUE(.pred$ndiff == .inner$ndiff))
  })

})
