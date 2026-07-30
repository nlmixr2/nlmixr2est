# Pareto k-hat: the tail-sensitive importance-sampling diagnostic.
#
# xi and the Kish effective sample size are both means over samples drawn FROM
# the proposal, so neither can see a tail the proposal does not visit.  k-hat
# estimates the tail index of the weight distribution and can.  These tests pin
# that the estimator actually recovers a known tail index, because a biased
# diagnostic is worse than no diagnostic.
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
                                                 isample = 300L, covMethod = "")))
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
    .m <- function() {
      ini({tka <- 0.45; tcl <- 1; tv <- 3.45; eta.cl ~ 0.1; add.sd <- 0.7})
      model({ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
             linCmt() ~ add(add.sd)})
    }
    .f <- suppressWarnings(nlmixr2(.m, nlmixr2data::theo_sd, "impmap",
                                   impmapControl(print = 0L, nIter = 6L,
                                                 isample = 300L, covMethod = "")))
    .E <- .f$env
    .bad <- which(.E$impPsisK > 0.7)
    # at least one subject is in the unreliable regime
    expect_gt(length(.bad), 0L)
    # ... and for those subjects the two in-sample statistics look FINE, which
    # is precisely why they could not be used to justify (or refute) the fix
    expect_true(all(abs(.E$impXi[.bad] - 1) < 0.15))
    expect_true(all((.E$impNeff / .E$impNsample)[.bad] > 0.9))
  })

  test_that("the k-hat alarm persists as isample grows (not a small-tail artifact)", {
    # A genuinely heavy tail is revealed MORE by extra samples; small-sample
    # noise would wash out.  Measured: max k-hat 2.49 (isample 300) -> 3.31
    # (1000) -> 3.76 (4000).
    skip_on_cran()
    .m <- function() {
      ini({tka <- 0.45; tcl <- 1; tv <- 3.45; eta.cl ~ 0.1; add.sd <- 0.7})
      model({ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
             linCmt() ~ add(add.sd)})
    }
    .maxK <- function(ns) {
      .f <- suppressWarnings(nlmixr2(.m, nlmixr2data::theo_sd, "impmap",
                                     impmapControl(print = 0L, nIter = 5L,
                                                   isample = ns, covMethod = "")))
      max(.f$env$impPsisK)
    }
    expect_gt(.maxK(300L), 1.0)
    expect_gt(.maxK(2000L), 1.0)
  })

})
