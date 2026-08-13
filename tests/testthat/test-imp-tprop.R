# Multivariate-t importance-sampling proposal (NONMEM DF).
#
# gamma can only make a Gaussian proposal WIDER; it cannot give it heavier
# tails.  When the target posterior's tails are heavier than the proposal's the
# importance weights have infinite variance -- detectable by Pareto k-hat and
# invisible to xi and the Kish ESS.  A t proposal has polynomial tails that
# dominate a Gaussian target's, which bounds the weights.
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

  # FIXTURE: THREE ETAs.  These tests need a Gaussian proposal that genuinely
  # FAILS, and that has to be measured, not assumed.  The previous fixture was a
  # single-ETA model, which reads k-hat about -1.4 for every subject and can never
  # fail -- the premise was unsatisfiable.  Measured max k-hat, seeds 42/43/44,
  # auto=FALSE, nIter=8, isample=300:
  #
  #   1 eta, full          -1.42 / -1.45 / -1.32    0 / 0 / 0 subjects above 0.7
  #   3 eta, full  df=0     0.84 /  0.74 /  0.66    2 / 1 / 0
  #   3 eta, full  df=20    0.65 /  0.40 /  0.30    0 / 0 / 0   clears every seed
  #
  # Sparse data (3 observations per subject) was tried and REJECTED: it fails on
  # every seed but no df clears it (df=3 still leaves 0.76), because there the
  # heavy tail is created by the DATA, and proposal shape cannot repair that --
  # the same reason the AUTO ladder withdraws an escalation that does not improve
  # k-hat.  A test asserting "the t proposal fixes it" needs a fixture where it
  # can.
  #
  # impSeed defaults to 42, so these fits are deterministic at the 0.84 row.
  .tModel <- function() {
    ini({tka <- 0.45; tcl <- 1; tv <- 3.45
         eta.ka ~ 0.1; eta.cl ~ 0.1; eta.v ~ 0.1; add.sd <- 0.7})
    model({ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
           linCmt() ~ add(add.sd)})
  }

  .fitDf <- function(df, nIter = 8L, isample = 300L) {
    suppressWarnings(nlmixr2(.tModel, nlmixr2data::theo_sd, "impmap",
                             impmapControl(print = 0L, nIter = nIter,
                                           isample = isample, covMethod = "",
                                           auto = FALSE, df = df, gammaRule = "floor")))
  }

  test_that("df control round-trips and defaults to the Gaussian proposal", {
    expect_equal(impmapControl()$df, 0)
    expect_equal(impmapControl(df = 5, gammaRule = "floor")$df, 5)
    expect_equal(do.call(impmapControl, impmapControl(df = 7, gammaRule = "floor"))$df, 7)
    expect_true("df" %in% .impmapIsControlNames)
  })

  test_that("a t proposal removes the tail failure a Gaussian one has", {
    # THE POINT OF THE FEATURE.  On plain theophylline the Gaussian proposal
    # leaves several subjects with infinite-variance weights (k-hat > 0.7); a t
    # proposal bounds them.  Measured max k-hat: 1.61 (Gaussian) vs -0.58
    # (df=20), -0.86 (df=5).
    .g <- .fitDf(0)
    .t <- .fitDf(20)
    expect_gt(max(.g$env$impPsisK), 0.7)      # premise: the Gaussian fails
    expect_gt(sum(.g$env$impPsisK > 0.7), 0L)
    expect_lt(max(.t$env$impPsisK), 0.7)      # and the t proposal fixes it
    expect_equal(sum(.t$env$impPsisK > 0.7), 0L)
  })

  test_that("the fix is cheap -- it barely costs effective sample size", {
    # The competing mechanism (widening a Gaussian via gamma) cost ESS
    # 0.95 -> 0.70 on this model.  df = 20 costs well under a percent, because
    # it changes the proposal's SHAPE rather than inflating its width.
    .g <- .fitDf(0)
    .t <- .fitDf(20)
    .essG <- mean(.g$env$impNeff / .g$env$impNsample)
    .essT <- mean(.t$env$impNeff / .t$env$impNsample)
    expect_gt(.essT, 0.95 * .essG)
  })

  test_that("df does not move the estimates, only the weights' behaviour", {
    # The importance weights divide by the proposal density, so changing the
    # proposal shape must not change what the EM converges to.
    .g <- .fitDf(0, nIter = 12L)
    .t <- .fitDf(5, nIter = 12L)
    expect_equal(unname(.t$theta), unname(.g$theta), tolerance = 0.02)
    expect_equal(.t$objf, .g$objf, tolerance = 0.5)
  })

  test_that("large df recovers the Gaussian normalizer", {
    # Ci gains (p/2)log(df/2) + lgamma(df/2) - lgamma((df+p)/2), which -> 0 as
    # df -> Inf.  If that correction were wrong (or omitted) the objective would
    # be offset by a constant, so agreement at large df is a direct check of the
    # normalizer rather than of the sampler.
    .g <- .fitDf(0, nIter = 5L, isample = 1000L)
    .t <- .fitDf(2000, nIter = 5L, isample = 1000L)
    expect_equal(.t$env$impObj, .g$env$impObj, tolerance = 0.15)
  })

  test_that("the t proposal composes with qrpem and with individual gamma", {
    .q <- suppressWarnings(nlmixr2(.tModel, nlmixr2data::theo_sd, "qrpem",
                                   qrpemControl(print = 0L, nIter = 6L,
                                                covMethod = "", df = 10)))
    expect_true(inherits(.q, "nlmixr2FitCore"))
    expect_true(all(is.finite(.q$env$impPsisK)))
    .i <- suppressWarnings(nlmixr2(.tModel, nlmixr2data::theo_sd, "impmap",
                                   impmapControl(print = 0L, nIter = 6L,
                                                 covMethod = "", df = 10,
                                                 gammaMethod = "individual", gammaRule = "floor")))
    expect_equal(.i$env$impGammaMethod, "individual")
    expect_true(all(is.finite(.i$env$impGammaInd)))
  })

  test_that("covMethod='imp' works under a t proposal", {
    # The covariance step builds its own proposal and must use the SAME df, or
    # the reweighted objective is not the one the fit converged on.
    .f <- suppressWarnings(nlmixr2(.tModel, nlmixr2data::theo_sd, "impmap",
                                   impmapControl(print = 0L, nIter = 8L,
                                                 covMethod = "imp", df = 10, gammaRule = "floor")))
    .cv <- as.matrix(.f$env$cov)
    expect_true(all(is.finite(.cv)))
    expect_true(all(diag(.cv) > 0))
  })

})

# Per-subject ISAMPLE (the NM7 Technical Guide's own remedy) vs the t proposal.
nmTest({


  # three ETAs, for the same reason .tModel has them: the tail-failure premise
  # below is unsatisfiable on a single-ETA model
  .isModel <- function() {
    ini({tka <- 0.45; tcl <- 1; tv <- 3.45
         eta.ka ~ 0.1; eta.cl ~ 0.1; eta.v ~ 0.1; add.sd <- 0.7})
    model({ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
           linCmt() ~ add(add.sd)})
  }

  test_that("isample accepts a per-subject vector and uses it", {
    expect_equal(impmapControl(isample = 300L, gammaRule = "floor")$isample, 300L)
    expect_equal(impmapControl(isample = c(100L, 200L), gammaRule = "floor")$isample, c(100L, 200L))
    expect_error(impmapControl(isample = 0L, gammaRule = "floor"))
    .d <- nlmixr2data::theo_sd
    .n <- length(unique(.d$ID))
    .iv <- rep(300L, .n); .iv[c(1L, 3L)] <- 900L
    # auto = FALSE: this tests that an explicit per-subject vector is PLUMBED
    # through, and AUTO's job is to reallocate the budget away from exactly such
    # a vector (measured: 836/912/905/... against the requested 900/300/300/...),
    # which would be testing AUTO rather than the plumbing.
    .f <- suppressWarnings(nlmixr2(.isModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 5L,
                                                 isample = .iv, covMethod = "",
                                                 auto = FALSE, gammaRule = "floor")))
    expect_equal(as.integer(.f$env$impNsampleInd), .iv)
  })

  test_that("a uniform isample vector matches the equivalent scalar", {
    .d <- nlmixr2data::theo_sd
    .n <- length(unique(.d$ID))
    .a <- suppressWarnings(nlmixr2(.isModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 5L,
                                                 isample = 300L, covMethod = "", gammaRule = "floor")))
    .b <- suppressWarnings(nlmixr2(.isModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 5L,
                                                 isample = rep(300L, .n),
                                                 covMethod = "", gammaRule = "floor")))
    expect_equal(.b$objf, .a$objf, tolerance = 1e-10)
    expect_equal(unname(.b$theta), unname(.a$theta), tolerance = 1e-10)
  })

  test_that("a heavier tail matches a 10x sample boost at a tenth the draws", {
    # ORIGINALLY this asserted the stronger claim that more samples CANNOT repair
    # a tail failure while a heavier tail can, on the numbers
    #
    #   Gaussian, uniform 300              max k-hat 1.611, 4 subjects bad
    #   Gaussian, bad subjects at 3000     max k-hat 3.276, 1 subject bad
    #   t proposal df=20, uniform 300      max k-hat -0.577, 0 bad
    #
    # That asymmetry NO LONGER REPRODUCES.  Re-measured (3 ETAs, auto=FALSE,
    # nIter=8, seed 42):
    #
    #            gaussian 300      boost bad->3000      df=20 at 300
    #   full     0.836 (2 bad)     0.678 (0 bad)        0.653 (0 bad)
    #   sparse   0.847 (6 bad)     0.859 (2 bad)        0.976 (2 bad)
    #
    # On full data BOTH remedies clear it; on sparse data NEITHER does, because
    # there the heavy tail is created by the DATA and no proposal repairs it.
    # So the sample-count-vs-shape asymmetry is not a property this fixture can
    # demonstrate any more.
    #
    # What survives, and is the practically useful half, is the COST: the t
    # proposal reaches the same max k-hat at the ORIGINAL sample count, i.e. for
    # a tenth of the draws spent on the boosted subjects.
    #
    # NOTE for whoever revisits AUTO: the retired claim is the stated rationale
    # for "sample count is deliberately not driven by Pareto k-hat"
    # (src/imp.cpp).  That rationale now rests on the sparse row alone.
    skip_on_cran()
    .d <- nlmixr2data::theo_sd             # 3-ETA fixture -- see .tModel above
    .n <- length(unique(.d$ID))
    .g <- suppressWarnings(nlmixr2(.isModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 8L,
                                                 isample = 300L, covMethod = "",
                                                 auto = FALSE, gammaRule = "floor")))
    .k0 <- .g$env$impPsisK
    expect_gt(sum(.k0 > 0.7), 0L)                     # premise
    .iv <- rep(300L, .n); .iv[.k0 > 0.7] <- 3000L
    .boost <- suppressWarnings(nlmixr2(.isModel, .d, "impmap",
                                       impmapControl(print = 0L, nIter = 8L,
                                                     isample = .iv, covMethod = "",
                                                     auto = FALSE, gammaRule = "floor")))
    .t <- suppressWarnings(nlmixr2(.isModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 8L,
                                                 isample = 300L, covMethod = "",
                                                 auto = FALSE, df = 20, gammaRule = "floor")))
    # the heavier tail is at least as good as 10x the draws ...
    expect_lte(max(.t$env$impPsisK), max(.boost$env$impPsisK) + 1e-8)
    # ... and clears the unreliable regime at the ORIGINAL sample count
    expect_lt(max(.t$env$impPsisK), 0.7)
    expect_equal(sum(.t$env$impPsisK > 0.7), 0L)
  })

})
