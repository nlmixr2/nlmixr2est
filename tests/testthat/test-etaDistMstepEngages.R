# The declared-distribution M-step (etaDistMstep=TRUE) is opt-in and every way
# it can decline to run is silent: no declaration in the model, a family the
# C++ dispatch does not implement, a copula block wider than a pair, or -- the
# one that actually happened -- metadata that cannot be built because
# rxEtaDistExpand() has already cleared the iniDf's etaDist column and replaced
# the copula block with rxCor.* thetas, so rxUiEtaDists() reports nothing on
# the expanded ui every estimator actually receives.
#
# With the M-step inert, etaDistMstep=TRUE and etaDistMstep=FALSE give
# BIT-IDENTICAL fits.  Estimates are therefore no evidence whatever that the
# mechanism engaged, and these tests assert engagement directly instead.
nmTest({
  test_that("the declaration survives rxEtaDistExpand() for the M-step", {
    .mod <- function() {
      ini({
        lclm <- 1.9; lv1m <- 1.8; tq <- 0.9; tv2 <- 4.2
        lclrv <- -2; lv1rv <- -2
        eta.cl + eta.v1 ~ c(1, 0.3, 1)
        dist(eta.cl) ~ dgamma(shape = 1 / exp(lclrv),
                              rate = 1 / (exp(lclrv) * exp(lclm)))
        dist(eta.v1) ~ dgamma(shape = 1 / exp(lv1rv),
                              rate = 1 / (exp(lv1rv) * exp(lv1m)))
        eta.q + eta.v2 ~ c(0.1, 0.01, 0.1)
        prop.sd <- 0.316
      })
      model({
        cl <- eta.cl; v <- eta.v1
        q <- exp(tq + eta.q); v2 <- exp(tv2 + eta.v2)
        linCmt() ~ prop(prop.sd)
      })
    }
    .u0 <- rxode2::rxUiDecompress(rxode2::rxode2(.mod))
    .d <- rxode2::rxUiEtaDists(.u0)
    expect_equal(nrow(.d), 2L)

    .stash <- nlmixr2est:::.etaDistDeclStash(.u0, .d)
    expect_false(is.null(.stash))
    expect_equal(.stash$name, c("eta.cl", "eta.v1"))
    # the copula pairs the second declared eta back to the first
    expect_equal(.stash$corWith, c(-1L, 0L))
    expect_equal(.stash$corTheta[2], "rxCor.eta.v1.eta.cl")

    # the expansion is lossy: this is exactly why the stash exists
    .u <- rxode2::rxUiDecompress(rxode2::rxEtaDistExpand(.u0))
    expect_equal(nrow(rxode2::rxUiEtaDists(.u)), 0L)
    expect_null(nlmixr2est:::.etaDistMstepCore(.u))

    # ... and with it, the metadata builds
    nlmixr2est:::.etaDistDeclSet(.u, .stash)
    .core <- nlmixr2est:::.etaDistMstepCore(.u)
    expect_false(is.null(.core))
    expect_equal(.core$fam, c(13L, 13L))          # dgamma
    expect_equal(.core$corWith, c(-1L, 0L))
    # rho comes from tanh() of the rxCor theta the expansion wrote, not from a
    # stale ini() value
    expect_equal(.core$rho[2], 0.3, tolerance = 1e-6)
    # gamma(shape = 1/exp(-2), rate = 1/(exp(-2)*exp(1.9)))
    expect_equal(.core$args[1, 1], exp(2), tolerance = 1e-6)
    expect_equal(.core$args[1, 2], 1 / (exp(-2) * exp(1.9)), tolerance = 1e-6)
  })

  test_that("the preProcess hook leaves the stash where the M-step finds it", {
    .mod <- function() {
      ini({
        lclm <- 1.9; tv <- 3.45; tka <- 0.45
        lclrv <- -2
        eta.cl ~ 1
        dist(eta.cl) ~ dgamma(shape = 1 / exp(lclrv),
                              rate = 1 / (exp(lclrv) * exp(lclm)))
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka); cl <- eta.cl; v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .u0 <- rxode2::rxUiDecompress(rxode2::rxode2(.mod))
    .r <- nlmixr2est:::.preProcessEtaDist(.u0, "saem", NULL,
                             saemControl(etaDistMstep = TRUE))
    expect_true(is.list(.r))
    .u <- rxode2::rxUiDecompress(.r$ui)
    # the hook must stash in `meta`: the ui environment, the control and an
    # extra iniDf column are all dropped before the estimator asks (the ui is
    # rebuilt, and the est method installs its own freshly-built control after
    # the hooks run).  Measured, not assumed.
    expect_false(is.null(nlmixr2est:::.etaDistDeclGet(.u)))
    expect_false(is.null(nlmixr2est:::.etaDistMstepCore(.u)))
  })

  test_that("saem's M-step actually runs, and moves the declared parameters", {
    .mod <- function() {
      ini({
        lclm <- 1.9; tv <- 3.45; tka <- 0.45
        lclrv <- -2
        eta.cl ~ 1
        dist(eta.cl) ~ dgamma(shape = 1 / exp(lclrv),
                              rate = 1 / (exp(lclrv) * exp(lclm)))
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka); cl <- eta.cl; v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .d <- nlmixr2data::theo_sd
    .ctl <- function(on) {
      # etaDistCorMstep defaults TRUE, so the copula closed form runs even when
      # the FAMILY M-step is off.  This test contrasts the family step, so the
      # "off" arm has to turn both off -- otherwise its counter is non-zero and
      # the contrast is not the one being asserted.  (That assertion only held
      # before because the copula step was itself inert: rxode2 bounded the
      # rxCor theta, .preProcessBoundedTransform renamed it to rxBoundedTr.*,
      # and the stash stopped resolving.)
      saemControl(nBurn = 5, nEm = 5, print = 0, covMethod = "",
                  etaDistMstep = on, etaDistCorMstep = on)
    }
    .fOn <- suppressWarnings(nlmixr2(.mod, .d, est = "saem", control = .ctl(TRUE)))
    # the counter is reset per fit, so this is THIS fit's count
    expect_gt(nlmixr2est:::saemEtaDistN_(), 0)

    .fOff <- suppressWarnings(nlmixr2(.mod, .d, est = "saem", control = .ctl(FALSE)))
    expect_equal(nlmixr2est:::saemEtaDistN_(), 0)

    # and it has to make a difference: with the M-step inert the two fits are
    # bit-identical, which is how this went unnoticed in the first place
    .on <- .fOn$parFixedDf[, "Estimate"]
    .off <- .fOff$parFixedDf[, "Estimate"]
    expect_false(isTRUE(all.equal(.on, .off)))
  })

  test_that("holding thetas out of the optimizer does not corrupt the parameter history", {
    # The history matrix is sized by the optimizer's npars while its column
    # NAMES come from nparsPrint, which also counts thetas that are reported but
    # held out of the optimizer.  When those disagree the recorded values shift
    # silently: with 2 of 5 thetas held out, the free thetas' values were written
    # under the held-out thetas' names and the last column was left an unwritten
    # NULL slot, which `[` then dropped.  muPrintActive() has to cover EVERY
    # hold-out for foceiOfvOptim to record nparsPrint columns instead of npars.
    #
    # Not hypothetical, and not specific to this feature: any future hold-out
    # that makes nparsPrint != npars hits the same path.
    .mod <- function() {
      ini({
        lclm <- log(3); lclrv <- -2; tka <- 0.45; tv <- 3.45
        eta.cl ~ 1
        dist(eta.cl) ~ dgamma(shape = 1 / exp(lclrv),
                              rate = 1 / (exp(lclrv) * exp(lclm)))
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka); cl <- eta.cl; v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .d <- nlmixr2data::theo_sd
    .fit <- function(on) {
      suppressWarnings(nlmixr2(.mod, .d, est = "focei",
        control = foceiControl(print = 1, etaDistMstep = on, covMethod = "")))
    }
    .on <- .fit(TRUE)
    .off <- .fit(FALSE)
    .rawOn <- get("parHistData", envir = .on$env)
    .rawOff <- get("parHistData", envir = .off$env)
    # every declared column actually carries values -- an unwritten slot comes
    # back NULL, which is what silently truncated the history
    expect_false(any(vapply(.rawOn, is.null, logical(1))))
    expect_equal(names(.rawOn), names(.rawOff))
    # and the processed history keeps every parameter, the held-out ones included
    expect_equal(colnames(.on$parHist), colnames(.off$parHist))
    expect_true(all(c("lclm", "lclrv", "add.sd") %in% colnames(.on$parHist)))
  })

  test_that("the copula correlation M-step is scale-invariant", {
    # This one line caused every saem divergence.  The old estimator was the
    # raw product-moment mean(z1*z2), which IS the constrained MLE for a
    # unit-variance bivariate normal -- but only when the draws actually have
    # unit variance.  They do not: the pooled latent spread starts around 2.15
    # and only settles near 0.8, so a true correlation of 0.3 came out as
    # 0.3 * mean(z^2) ~ 1.38 and clamped to 0.999 on the FIRST M-step.  It then
    # self-sustained (0.999 * mean(z^2) is itself ~0.999), and a clamped rho
    # makes the copula partner's latent w_k = rho*z_j + sqrt(1-rho^2)*z_k
    # numerically equal to z_j -- so BOTH declared families end up fitted to the
    # same draws.
    #
    # The sample correlation is bounded by construction and identical to the
    # product-moment when the draws are standard, so these assertions pin the
    # property that actually matters: the answer must not move when the draws
    # are scaled.
    set.seed(42)
    .n <- 4000
    .rho <- 0.3
    .z1 <- stats::rnorm(.n)
    .z2 <- .rho * .z1 + sqrt(1 - .rho^2) * stats::rnorm(.n)

    # standard draws: recovers rho
    expect_equal(nlmixr2est:::rxEtaDistCorTest_(.z1, .z2), .rho, tolerance = 0.05)

    # the case that broke it -- an unmixed chain, spread ~2
    expect_equal(nlmixr2est:::rxEtaDistCorTest_(2 * .z1, 2 * .z2), .rho, tolerance = 0.05)
    expect_lt(nlmixr2est:::rxEtaDistCorTest_(2 * .z1, 2 * .z2), 0.9)   # never near the clamp

    # and the shrunk case, spread ~0.8, which is where a settled fit sits
    expect_equal(nlmixr2est:::rxEtaDistCorTest_(0.8 * .z1, 0.8 * .z2), .rho, tolerance = 0.05)

    # scale invariance outright, including unequal scales
    expect_equal(nlmixr2est:::rxEtaDistCorTest_(3 * .z1, 0.5 * .z2),
                 nlmixr2est:::rxEtaDistCorTest_(.z1, .z2), tolerance = 1e-8)

    # a location shift must not matter either
    expect_equal(nlmixr2est:::rxEtaDistCorTest_(.z1 + 5, .z2 - 2),
                 nlmixr2est:::rxEtaDistCorTest_(.z1, .z2), tolerance = 1e-8)

    # still bounded for a genuinely near-perfect correlation
    expect_lte(nlmixr2est:::rxEtaDistCorTest_(.z1, .z1), 0.999)
  })
  test_that("the M-step actually FIRES in a saem fit, not just maps", {
    # The tests above check the R-side core builds from a stashed declaration.
    # None of them ran a fit, and that gap hid a real regression: rxode2 gave
    # the copula theta finite bounds, .preProcessBoundedTransform wrapped it and
    # renamed it to rxBoundedTr.rxCor.*, the stash's recorded name stopped
    # resolving, and the M-step went SILENTLY inert -- etaDistMstep=TRUE
    # produced estimates identical to FALSE, with only a runInfo warning to say
    # so.  Assert engagement against the counter, in a real fit.
    skip_if_not(file.exists("~/src/gamma_indpar/gamma_clv1.dat"))
    .d <- utils::read.table("~/src/gamma_indpar/gamma_clv1.dat", skip = 1,
                            header = TRUE)
    names(.d) <- c("ID", "TIME", "AMT", "RATE", "EVID", "MDV", "DV", "IPRED")
    .d$CMT <- 1
    .d <- .d[, c("ID", "TIME", "AMT", "RATE", "EVID", "DV", "CMT")]
    .m <- function() {
      ini({
        lclm <- 1.9; lv1m <- 1.8; tq <- 0.9; tv2 <- 4.2
        lclrv <- -2; lv1rv <- -2
        eta.cl + eta.v1 ~ c(1, 0.3, 1)
        dist(eta.cl) ~ dgamma(shape = 1 / exp(lclrv),
                              rate = 1 / (exp(lclrv) * exp(lclm)))
        dist(eta.v1) ~ dgamma(shape = 1 / exp(lv1rv),
                              rate = 1 / (exp(lv1rv) * exp(lv1m)))
        eta.q + eta.v2 ~ c(0.1, 0.01, 0.1)
        prop.sd <- 0.316
      })
      model({
        cl <- eta.cl; v <- eta.v1
        q <- exp(tq + eta.q); v2 <- exp(tv2 + eta.v2)
        linCmt() ~ prop(prop.sd)
      })
    }
    .f <- suppressWarnings(nlmixr2(.m, .d, est = "saem",
      control = saemControl(nBurn = 25, nEm = 25, print = 0L, covMethod = "",
                            seed = 99, calcTables = FALSE,
                            etaDistMstep = TRUE)))
    expect_gt(nlmixr2est:::saemEtaDistN_(), 0L)
    expect_false(any(grepl("etaDistMstep=TRUE was requested",
                           as.character(.f$runInfo), fixed = TRUE)))
  })

})
