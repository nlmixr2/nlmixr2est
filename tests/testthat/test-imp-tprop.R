# Multivariate-t importance-sampling proposal (NONMEM DF).
#
# gamma can only make a Gaussian proposal WIDER; it cannot give it heavier
# tails.  When the target posterior's tails are heavier than the proposal's the
# importance weights have infinite variance -- detectable by Pareto k-hat and
# invisible to xi and the Kish ESS.  A t proposal has polynomial tails that
# dominate a Gaussian target's, which bounds the weights.
nmTest({

  .tModel <- function() {
    ini({tka <- 0.45; tcl <- 1; tv <- 3.45; eta.cl ~ 0.1; add.sd <- 0.7})
    model({ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
           linCmt() ~ add(add.sd)})
  }

  .fitDf <- function(df, nIter = 8L, isample = 300L) {
    suppressWarnings(nlmixr2(.tModel, nlmixr2data::theo_sd, "impmap",
                             impmapControl(print = 0L, nIter = nIter,
                                           isample = isample, covMethod = "",
                                           df = df)))
  }

  test_that("df control round-trips and defaults to the Gaussian proposal", {
    expect_equal(impmapControl()$df, 0)
    expect_equal(impmapControl(df = 5)$df, 5)
    expect_equal(do.call(impmapControl, impmapControl(df = 7))$df, 7)
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
                                                 gammaMethod = "individual")))
    expect_equal(.i$env$impGammaMethod, "individual")
    expect_true(all(is.finite(.i$env$impGammaInd)))
  })

  test_that("covMethod='imp' works under a t proposal", {
    # The covariance step builds its own proposal and must use the SAME df, or
    # the reweighted objective is not the one the fit converged on.
    .f <- suppressWarnings(nlmixr2(.tModel, nlmixr2data::theo_sd, "impmap",
                                   impmapControl(print = 0L, nIter = 8L,
                                                 covMethod = "imp", df = 10)))
    .cv <- as.matrix(.f$env$cov)
    expect_true(all(is.finite(.cv)))
    expect_true(all(diag(.cv) > 0))
  })

})

# Per-subject ISAMPLE (the NM7 Technical Guide's own remedy) vs the t proposal.
nmTest({

  .isModel <- function() {
    ini({tka <- 0.45; tcl <- 1; tv <- 3.45; eta.cl ~ 0.1; add.sd <- 0.7})
    model({ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
           linCmt() ~ add(add.sd)})
  }

  test_that("isample accepts a per-subject vector and uses it", {
    expect_equal(impmapControl(isample = 300L)$isample, 300L)
    expect_equal(impmapControl(isample = c(100L, 200L))$isample, c(100L, 200L))
    expect_error(impmapControl(isample = 0L))
    .d <- nlmixr2data::theo_sd
    .n <- length(unique(.d$ID))
    .iv <- rep(300L, .n); .iv[c(1L, 3L)] <- 900L
    .f <- suppressWarnings(nlmixr2(.isModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 5L,
                                                 isample = .iv, covMethod = "")))
    expect_equal(as.integer(.f$env$impNsampleInd), .iv)
  })

  test_that("a uniform isample vector matches the equivalent scalar", {
    .d <- nlmixr2data::theo_sd
    .n <- length(unique(.d$ID))
    .a <- suppressWarnings(nlmixr2(.isModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 5L,
                                                 isample = 300L, covMethod = "")))
    .b <- suppressWarnings(nlmixr2(.isModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 5L,
                                                 isample = rep(300L, .n),
                                                 covMethod = "")))
    expect_equal(.b$objf, .a$objf, tolerance = 1e-10)
    expect_equal(unname(.b$theta), unname(.a$theta), tolerance = 1e-10)
  })

  test_that("more samples does NOT repair tail failure -- only a heavier tail does", {
    # The guide's remedy treats the symptom; the tutorial's treats the cause.
    # Infinite weight variance is a property of the proposal/target TAIL RATIO,
    # so extra draws from a too-light proposal simply reach further into the
    # region that breaks the estimator and can make k-hat WORSE.  Measured:
    #
    #   Gaussian, uniform 300              max k-hat 1.611, 4 subjects bad
    #   Gaussian, bad subjects at 3000     max k-hat 3.276, 1 subject bad
    #   t proposal df=20, uniform 300      max k-hat -0.577, 0 bad
    skip_on_cran()
    .d <- nlmixr2data::theo_sd
    .n <- length(unique(.d$ID))
    .g <- suppressWarnings(nlmixr2(.isModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 8L,
                                                 isample = 300L, covMethod = "")))
    .k0 <- .g$env$impPsisK
    expect_gt(sum(.k0 > 0.7), 0L)                     # premise
    .iv <- rep(300L, .n); .iv[.k0 > 0.7] <- 3000L
    .boost <- suppressWarnings(nlmixr2(.isModel, .d, "impmap",
                                       impmapControl(print = 0L, nIter = 8L,
                                                     isample = .iv, covMethod = "")))
    .t <- suppressWarnings(nlmixr2(.isModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 8L,
                                                 isample = 300L, covMethod = "",
                                                 df = 20)))
    # 10x the samples does not clear the unreliable regime ...
    expect_gt(max(.boost$env$impPsisK), 0.7)
    # ... while a heavier-tailed proposal does, at the ORIGINAL sample count
    expect_lt(max(.t$env$impPsisK), 0.7)
  })

})
