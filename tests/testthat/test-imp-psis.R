# Pareto k-hat: the tail-sensitive importance-sampling diagnostic.
#
# xi and the Kish effective sample size are both means over samples drawn FROM
# the proposal, so neither can see a tail the proposal does not visit.  k-hat
# estimates the tail index of the weight distribution and can.  These tests pin
# that the estimator actually recovers a known tail index, because a biased
# diagnostic is worse than no diagnostic.
# NOTE ON gammaRule.  These tests exercise the TAIL machinery -- the t proposal
# (df), Pareto k-hat, and AUTO's df escalation -- all of which need a Gaussian
# proposal that actually FAILS in order to have anything to repair.  The default
# rule is now "target", which drives xi onto iaccept and in doing so repairs the
# tail itself: on the 3-ETA theophylline fixture it takes max k-hat from 0.836 to
# about -0.4, leaving these premises unsatisfiable.  So they pin
# gammaRule = "floor" deliberately.
#
# That overlap is a real consequence of the default change, not a test artifact:
# with "target" as the default the df/AUTO tail machinery is a secondary safety
# net rather than the primary remedy.
nmTest({

  .rgpd <- function(n, k) ((1 - stats::runif(n))^(-k) - 1) / k

  test_that("k-hat recovers a known generalized-Pareto tail index", {
    .est <- function(k, seeds = 1:5) {
      vapply(seeds, function(s) {
        set.seed(100 + s)
        .impPsisK(.rgpd(4000, k))
      }, numeric(1))
    }
    for (.k in c(0.2, 0.5, 0.8, 1.2)) {
      .m <- mean(.est(.k))
      expect_true(is.finite(.m))
      # within 0.15 of truth: comfortably tighter than the 0.5/0.7 decision
      # thresholds the statistic is actually used against
      expect_lt(abs(.m - .k), 0.15)
    }
  })

  test_that("k-hat separates light from heavy tails", {
    set.seed(11)
    # exponential weights: tail index 0
    expect_lt(.impPsisK(stats::rexp(4000)), 0.3)
    # heavy: k = 1 has infinite variance AND infinite mean
    expect_gt(mean(vapply(1:3, function(s) {
      set.seed(200 + s); .impPsisK(.rgpd(4000, 1.0))
    }, numeric(1))), 0.7)
  })

  test_that("k-hat is scale invariant and degrades gracefully", {
    set.seed(5)
    .w <- .rgpd(2000, 0.6)
    # scale invariance matters: the stored weights carry a per-subject mixture
    # responsibility factor that must not shift the diagnostic
    expect_equal(.impPsisK(.w), .impPsisK(.w * 1e6), tolerance = 1e-8)
    expect_equal(.impPsisK(.w), .impPsisK(.w * 1e-6), tolerance = 1e-8)
    # too few weights to fit a tail -> NA rather than a wrong number
    expect_true(is.na(.impPsisK(.w[1:10])))
    expect_true(is.na(.impPsisK(numeric(0))))
    # non-finite and non-positive weights are dropped, not propagated
    expect_true(is.finite(.impPsisK(c(.w, NA, Inf, 0, -1))))
  })

  test_that("k-hat agrees with the loo reference implementation", {
    skip_if_not_installed("loo")
    for (.k in c(0.3, 0.7, 1.1)) {
      set.seed(321)
      .w <- .rgpd(4000, .k)
      .mine <- .impPsisK(.w)
      .ref <- suppressWarnings(
        loo::psis(matrix(log(.w), ncol = 1), r_eff = NA)$diagnostics$pareto_k)
      expect_lt(abs(.mine - .ref), 0.1)
    }
  })

  test_that("a fit exposes per-subject k-hat alongside xi and the Kish ESS", {
    .m <- function() {
      ini({tka <- 0.45; tcl <- 1; tv <- 3.45; eta.cl ~ 0.1; add.sd <- 0.7})
      model({ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
             linCmt() ~ add(add.sd)})
    }
    .d <- nlmixr2data::theo_sd
    .f <- suppressWarnings(nlmixr2(.m, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 6L,
                                                 isample = 300L, covMethod = "", gammaRule = "floor")))
    .E <- .f$env
    expect_equal(length(.E$impPsisK), length(unique(.d$ID)))
    # all three diagnostics are present and per-subject
    expect_equal(length(.E$impXi), length(.E$impPsisK))
    expect_equal(length(.E$impNeff), length(.E$impPsisK))
    expect_true(all(is.finite(.E$impPsisK)))
  })

  test_that("k-hat sees tail failure that xi and the Kish ESS cannot", {
    # THE MOTIVATING MEASUREMENT.  On plain theophylline with a Gaussian
    # (Laplace) proposal, some subjects have importance weights with infinite
    # variance -- k-hat well above the 0.7 reliability threshold -- while xi
    # sits at ~1.0 and the Kish effective-sample fraction at ~0.99 for those
    # SAME subjects.  Neither in-sample statistic can see it, because the
    # offending mass is in a tail the proposal rarely visits and no single
    # drawn weight dominates.
    #
    # This is the premise for the t-distribution proposal (NONMEM's DF): the
    # fix has to make the proposal's tails HEAVIER than the target's, which
    # widening a Gaussian by gamma cannot do.
    # THREE etas, and auto = FALSE.  One eta on theophylline has no tail failure
    # left to see (max k-hat about -1.8), and `auto` is on by default -- it would
    # escalate the proposal and repair the very thing this test exists to
    # observe.  Both matter: the fixture has to fail, and the failure has to be
    # left alone.
    .m <- function() {
      ini({tka <- 0.45; tcl <- 1; tv <- 3.45
           eta.ka ~ 0.3; eta.cl ~ 0.1; eta.v ~ 0.1; add.sd <- 0.7})
      model({ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
             linCmt() ~ add(add.sd)})
    }
    .f <- suppressWarnings(nlmixr2(.m, nlmixr2data::theo_sd, "impmap",
                                   impmapControl(print = 0L, nIter = 6L,
                                                 isample = 300L, covMethod = "",
                                                 auto = FALSE, gammaRule = "floor")))
    .E <- .f$env
    .bad <- which(.E$impPsisK > 0.7)
    # at least one subject is in the unreliable regime
    expect_gt(length(.bad), 0L)
    # ... and for those subjects the two in-sample statistics look FINE, which
    # is precisely why they could not be used to justify (or refute) the fix
    expect_true(all(abs(.E$impXi[.bad] - 1) < 0.15))
    expect_true(all((.E$impNeff / .E$impNsample)[.bad] > 0.9))
  })

  test_that("a STRUCTURAL tail failure persists as isample grows; noise washes out", {
    # The original form of this test asserted max k-hat > 1 at isample 300 AND
    # 2000 on one-eta theophylline, citing 2.49 -> 3.31 -> 3.76 as k-hat GROWING
    # with more draws.  That growth was itself a symptom of the pooling bug: a
    # well-behaved sampler improves with more draws, it does not degrade.  The
    # claim worth testing now is the discriminating one -- a tail failure the
    # DATA creates survives more sampling, while one that is sampling noise does
    # not.
    skip_on_cran()
    .m <- function() {
      ini({tka <- 0.45; tcl <- 1; tv <- 3.45
           eta.ka ~ 0.3; eta.cl ~ 0.1; eta.v ~ 0.1; add.sd <- 0.7})
      model({ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
             linCmt() ~ add(add.sd)})
    }
    # nobs < neta: two observations against three etas, so the individual
    # posterior is not identified and the heavy tail is a property of the model,
    # not of the proposal.
    set.seed(42)
    .sparse <- do.call(rbind, lapply(split(nlmixr2data::theo_sd,
                                           nlmixr2data::theo_sd$ID), function(d) {
      .obs <- d[d$EVID == 0, , drop = FALSE]
      rbind(d[d$EVID != 0, , drop = FALSE],
            .obs[sort(sample(seq_len(nrow(.obs)), 2L)), , drop = FALSE])
    }))
    .maxK <- function(dat, ns) {
      .f <- suppressWarnings(nlmixr2(.m, dat, "impmap",
                                     impmapControl(print = 0L, nIter = 5L,
                                                   isample = ns, covMethod = "",
                                                   auto = FALSE, gammaRule = "floor")))
      max(.f$env$impPsisK)
    }
    # structural: still unreliable at 300 and at 2000 (and the isample = 8000
    # reference for this fixture reads 0.794, i.e. it never clears)
    expect_gt(.maxK(.sparse, 300L), 0.7)
    expect_gt(.maxK(.sparse, 2000L), 0.7)
    # and the contrast: on data that DOES identify each subject, an apparent
    # failure at 300 washes out by 2000 -- which is why "k-hat > 0.7 once" is
    # not on its own evidence of a structural problem
    expect_gt(.maxK(nlmixr2data::theo_sd, 300L), 0.7)
    expect_lt(.maxK(nlmixr2data::theo_sd, 2000L), 0.7)
  })

})
