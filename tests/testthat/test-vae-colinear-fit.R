## End-to-end colinearity handling in the VAE covariate M-step: the hysteresis
## and near-tie passes only exist inside a real training run, so these need
## actual fits.  Weekly, not essential (see .slowBatches in tests/testthat.R);
## the pure helpers are covered by test-vae-colinear.R.

nmTest({

  ## CL genuinely depends on WT; LBM is a ~0.99 correlated decoy that carries no
  ## independent signal.  This is the situation the feature exists for: the two
  ## score within noise of one another, so the winner is near-arbitrary.
  .colinData <- function(nid = 40L, seed = 9L) {
    .testSeed(seed)
    wt <- round(stats::runif(nid, 50, 100), 1)
    lbm <- round(0.75 * wt + stats::rnorm(nid, 0, 1.2), 1)
    ctr <- stats::median(wt)
    cl <- 2.7 * exp(0.75 * log(wt / ctr) + stats::rnorm(nid, 0, 0.15))
    v <- 31 * exp(stats::rnorm(nid, 0, 0.1))
    ka <- 1.5 * exp(stats::rnorm(nid, 0, 0.3))
    tms <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 24)
    do.call(rbind, lapply(seq_len(nid), function(i) {
      ke <- cl[i] / v[i]
      f <- 320 / v[i] * ka[i] / (ka[i] - ke) * (exp(-ke * tms) - exp(-ka[i] * tms))
      rbind(data.frame(ID = i, TIME = 0, AMT = 320, EVID = 1, DV = 0,
                       WT = wt[i], LBM = lbm[i]),
            data.frame(ID = i, TIME = tms, AMT = 0, EVID = 0,
                       DV = f + stats::rnorm(length(tms), 0, 0.25),
                       WT = wt[i], LBM = lbm[i]))
    }))
  }

  .colinModel <- function() {
    ini({
      tka <- log(1.5); tcl <- log(2.7); tv <- log(31)
      eta.ka ~ 0.3; eta.cl ~ 0.2; eta.v ~ 0.1
      add.sd <- 0.3
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  .ctl <- function(...) {
    vaeControl(iters = 60L, itersBurnIn = 15L, calcTables = FALSE, ...)
  }

  .fitQuiet <- function(d, ...) {
    suppressMessages(suppressWarnings(
      nlmixr2(.colinModel, d, est = "vae", control = .ctl(...))))
  }

  test_that("the colinear pass runs and reports near-tied alternates", {
    skip_on_cran()
    d <- .colinData()
    ## the design really is clustered, so a counter reaching zero would mean the
    ## mechanism never engaged rather than that there was nothing to do
    cv <- vaeCovariates(d, warn = FALSE)
    expect_true(.vaeClusterBinds(cv$cluster, cv$group))

    ## a wide window so the report is exercised deterministically; the default is
    ## much tighter and correctly stays quiet on this fit
    f <- .fitQuiet(d, covSelectAltTol = 20)
    cd <- f$vae$colinear

    ## mechanism ran
    expect_gt(cd$nColinTest, 0L)
    ## and produced the diagnostic
    alt <- cd$alternates
    expect_s3_class(alt, "data.frame")
    expect_gt(nrow(alt), 0L)
    expect_identical(names(alt), c("param", "covariate", "alternate", "delta"))
    ## every reported alternate is a DIFFERENT covariate.  Alternate shapes of
    ## one covariate live in the same cluster and are near-tied on almost every
    ## model; admitting them would make this fire universally.
    expect_false(any(sub("_.*", "", alt$covariate) == sub("_.*", "", alt$alternate)))
    ## the true covariate is what got selected, and the decoy is the alternative
    expect_true(all(grepl("^WT", alt$covariate)))
    expect_true(all(grepl("^LBM", alt$alternate)))
    ## and the user is told, rather than having to go looking
    expect_match(f$runInfo, "near-tied colinear", all = FALSE)
  })

  test_that("the default window stays quiet on a well-determined fit", {
    skip_on_cran()
    f <- .fitQuiet(.colinData())
    ## a 0.99 decoy still fits measurably worse when the effect is strong, so at
    ## the default tolerance there is nothing to report and no note is raised
    expect_identical(nrow(f$vae$colinear$alternates), 0L)
    expect_false(any(grepl("near-tied colinear", f$runInfo)))
    ## the pass still ran -- silence here means "nothing was close", not "the
    ## mechanism was never reached"
    expect_gt(f$vae$colinear$nColinTest, 0L)
  })

  test_that("covSelectColinear=FALSE makes the mechanism unreachable", {
    skip_on_cran()
    d <- .colinData()
    fOff <- .fitQuiet(d, covSelectColinear = FALSE, covSelectAltTol = 20)
    ## not merely "reported nothing" -- never scored a single swap
    expect_identical(fOff$vae$colinear$nColinTest, 0L)
    expect_identical(fOff$vae$colinear$nColinHold, 0L)
    expect_identical(nrow(fOff$vae$colinear$alternates), 0L)
    expect_false(any(grepl("near-tied colinear", fOff$runInfo)))
  })

  ## A deliberately UNSTABLE fixture: a small study, a modest covariate effect
  ## and a ~0.999 decoy, fitted against the unsmoothed posterior means.  The
  ## winner genuinely flips between iterations here, which is the only situation
  ## in which the hold can fire at all.
  .flipData <- function(nid = 30L, seed = 7L) {
    .testSeed(seed)
    wt <- round(stats::runif(nid, 50, 100), 1)
    lbm <- round(0.75 * wt + stats::rnorm(nid, 0, 0.3), 1)
    ctr <- stats::median(wt)
    cl <- 2.7 * exp(0.5 * log(wt / ctr) + stats::rnorm(nid, 0, 0.2))
    v <- 31 * exp(stats::rnorm(nid, 0, 0.12))
    ka <- 1.5 * exp(stats::rnorm(nid, 0, 0.35))
    tms <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 24)
    do.call(rbind, lapply(seq_len(nid), function(i) {
      ke <- cl[i] / v[i]
      f <- 320 / v[i] * ka[i] / (ka[i] - ke) * (exp(-ke * tms) - exp(-ka[i] * tms))
      rbind(data.frame(ID = i, TIME = 0, AMT = 320, EVID = 1, DV = 0,
                       WT = wt[i], LBM = lbm[i]),
            data.frame(ID = i, TIME = tms, AMT = 0, EVID = 0,
                       DV = f + stats::rnorm(length(tms), 0, 0.4),
                       WT = wt[i], LBM = lbm[i]))
    }))
  }

  .flipModel <- function() {
    ini({
      tka <- log(1.5); tcl <- log(2.7); tv <- log(31)
      eta.ka ~ 0.35; eta.cl ~ 0.25; eta.v ~ 0.12
      add.sd <- 0.4
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("hysteresis actually holds an incumbent, and only because of the margin", {
    skip_on_cran()
    d <- .flipData()
    .fit <- function(h) {
      suppressMessages(suppressWarnings(
        nlmixr2(.flipModel, d, est = "vae",
                control = vaeControl(iters = 90L, itersBurnIn = 15L,
                                     calcTables = FALSE,
                                     covSelectSmooth = FALSE,
                                     covSelectHysteresis = h))))
    }
    wide <- .fit(20)
    ## the mechanism did not merely run -- it retained an incumbent the search
    ## had displaced
    expect_gt(wide$vae$colinear$nColinTest, 0L)
    expect_gt(wide$vae$colinear$nColinHold, 0L)
    ## and the MARGIN is why.  With no margin the same fixture holds nothing,
    ## so an implementation that always reverted could not pass both halves.
    none <- .fit(0)
    expect_gt(none$vae$colinear$nColinTest, 0L)
    expect_identical(none$vae$colinear$nColinHold, 0L)
  })

  ## eta.cl and eta.v are correlated BY CONSTRUCTION, which is what the
  ## cross-parameter refinement needs to have anything to arbitrate
  .phiData <- function(nid = 45L, seed = 11L) {
    .testSeed(seed)
    wt <- round(stats::runif(nid, 50, 100), 1)
    ctr <- stats::median(wt)
    ## correlated etas via a Cholesky factor rather than MASS::mvrnorm, which
    ## would add a package dependency for one fixture
    .sig <- matrix(c(0.09, 0.075, 0.075, 0.09), 2, 2)
    z <- matrix(stats::rnorm(nid * 2L), nid, 2L) %*% chol(.sig)
    cl <- 2.7 * exp(0.7 * log(wt / ctr) + z[, 1])
    v <- 31 * exp(z[, 2])
    ka <- 1.5 * exp(stats::rnorm(nid, 0, 0.3))
    tms <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 24)
    do.call(rbind, lapply(seq_len(nid), function(i) {
      ke <- cl[i] / v[i]
      f <- 320 / v[i] * ka[i] / (ka[i] - ke) * (exp(-ke * tms) - exp(-ka[i] * tms))
      rbind(data.frame(ID = i, TIME = 0, AMT = 320, EVID = 1, DV = 0, WT = wt[i]),
            data.frame(ID = i, TIME = tms, AMT = 0, EVID = 0,
                       DV = f + stats::rnorm(length(tms), 0, 0.25), WT = wt[i]))
    }))
  }

  .phiBlock <- function() {
    ini({
      tka <- log(1.5); tcl <- log(2.7); tv <- log(31)
      eta.ka ~ 0.3
      eta.cl + eta.v ~ c(0.09, 0.07, 0.09)
      add.sd <- 0.3
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  .phiDiag <- function() {
    ini({
      tka <- log(1.5); tcl <- log(2.7); tv <- log(31)
      eta.ka ~ 0.3; eta.cl ~ 0.09; eta.v ~ 0.09
      add.sd <- 0.3
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  .phiCtl <- function(...) {
    vaeControl(iters = 60L, itersBurnIn = 15L, calcTables = FALSE,
               covSelectPhiJoin = 0.5, covSelectPhiLeave = 0.4, ...)
  }

  test_that("a correlated omega lets the cross-parameter refinement run", {
    skip_on_cran()
    f <- suppressMessages(suppressWarnings(
      nlmixr2(.phiBlock, .phiData(), est = "vae", control = .phiCtl())))
    cd <- f$vae$colinear
    ## the gate opened, because the model declares the correlation
    expect_true(cd$omOff)
    ## groups formed and moves were evaluated -- "the mechanism ran"
    expect_gt(cd$nPhiPair, 0L)
    expect_gt(cd$nPhiTest, 0L)
    ## no advisory: the refinement was not skipped, so there is nothing to advise
    expect_false(any(grepl("declare an omega block", f$runInfo)))
  })

  test_that("a diagonal omega detects the groups but skips the refinement", {
    skip_on_cran()
    f <- suppressMessages(suppressWarnings(
      nlmixr2(.phiDiag, .phiData(), est = "vae", control = .phiCtl())))
    cd <- f$vae$colinear
    expect_false(cd$omOff)
    ## the dims ARE correlated and that is reported...
    expect_gt(cd$nPhiPair, 0L)
    ## ...but with a diagonal omega the objective is separable across dims, so
    ## not one move is scored -- this is a proof of the gate, not of an empty
    ## result
    expect_identical(cd$nPhiTest, 0L)
    expect_identical(cd$nPhiMove, 0L)
    ## and the modeler is told what would enable it
    expect_match(f$runInfo, "declare an omega block", all = FALSE)
  })

  test_that("covSelectColinear=FALSE reaches the phi refinement too", {
    skip_on_cran()
    f <- suppressMessages(suppressWarnings(
      nlmixr2(.phiBlock, .phiData(), est = "vae",
              control = vaeControl(iters = 60L, itersBurnIn = 15L,
                                   calcTables = FALSE,
                                   covSelectColinear = FALSE))))
    cd <- f$vae$colinear
    ## not even the grouping is computed
    expect_identical(cd$nPhiPair, 0L)
    expect_identical(cd$nPhiTest, 0L)
    expect_false(any(grepl("declare an omega block", f$runInfo)))
  })

  test_that("every covSelectPhiCor source drives the grouping in a real fit", {
    skip_on_cran()
    ## Each source builds the correlation matrix from a different quantity in
    ## C++ -- the smoothed sufficient statistic, the raw posterior means, or the
    ## means less the fitted covariate centers.  Only the default runs unless
    ## this asks for the others, so the two alternative branches would otherwise
    ## never execute.
    d <- .phiData(seed = 31L)
    for (src in c("suffStat", "mu", "resid")) {
      f <- suppressMessages(suppressWarnings(
        nlmixr2(.phiBlock, d, est = "vae",
                control = .phiCtl(covSelectPhiCor = src))))
      ## the branch ran and produced a grouping
      expect_gt(f$vae$colinear$nPhiPair, 0L, label = src)
      expect_gt(f$vae$colinear$nPhiTest, 0L, label = src)
    }
  })

  test_that("a wider sticky band never joins fewer pairs", {
    skip_on_cran()
    ## covSelectPhiLeave below covSelectPhiJoin is what keeps a pair joined when
    ## its correlation dips into the gap.  Widening the band can therefore only
    ## ever ADD joined pair-iterations, never remove them, whatever the
    ## correlation trajectory happens to be -- which is the part of stickiness
    ## that can be asserted without depending on that trajectory.
    ##
    ## KNOWN GAP: the specific case "correlation crossed below join and the pair
    ## stayed" is not directly asserted.  It needs a fixture whose correlation
    ## dips into the band and stays there, and on this fixture the effect is a
    ## single iteration (13 joined pair-iterations against 12).  A test hinging
    ## on that margin would be decided by differences at the 1e-14 level, which
    ## this fit already shows across thread counts, so it would be flaky rather
    ## than informative.  $vae$colinear$phiPairOn exposes the adjacency for
    ## anyone constructing a sharper fixture later.
    d <- .phiData(seed = 31L)
    .pairs <- function(lv) {
      f <- suppressMessages(suppressWarnings(
        nlmixr2(.phiBlock, d, est = "vae",
                control = vaeControl(iters = 80L, itersBurnIn = 15L,
                                     calcTables = FALSE,
                                     covSelectPhiCor = "mu",
                                     covSelectPhiJoin = 0.78,
                                     covSelectPhiLeave = lv))))
      list(n = f$vae$colinear$nPhiPair, on = f$vae$colinear$phiPairOn)
    }
    none <- .pairs(0.78)     # no band: leave == join
    wide <- .pairs(0.50)
    expect_gte(wide$n, none$n)
    ## the adjacency is surfaced, square, and symmetric -- a pair is a pair
    ## whichever way round it is read
    expect_true(is.matrix(wide$on))
    expect_identical(nrow(wide$on), ncol(wide$on))
    expect_identical(wide$on, t(wide$on))
  })

  test_that("a bounded intercept is held, and an emptied support survives it", {
    skip_on_cran()
    ## When a joint candidate would push a population intercept past its ini()
    ## bound, that intercept is HELD at the bound and the candidate re-scored --
    ## the per-dim path clamps AFTER scoring, which is harmless there but not
    ## here, where the score decides between moves.
    ##
    ## This also covers the case that used to break: holding the intercept drops
    ## the intercept column, and a DROP move can leave a support with no columns
    ## at all.  Nothing else in this file declares bounds, so without this test
    ## the held-intercept path never runs.
    .bounded <- function() {
      ini({
        tka <- log(1.5)
        tcl <- c(-Inf, log(2.7), log(2.0))   # upper bound below the truth
        tv <- c(-Inf, log(31), log(25))      # so both intercepts clamp
        eta.ka ~ 0.3
        eta.cl + eta.v ~ c(0.09,
                           0.07, 0.09)
        add.sd <- 0.3
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    f <- suppressMessages(suppressWarnings(
      nlmixr2(.bounded, .phiData(seed = 31L), est = "vae", control = .phiCtl())))
    cd <- f$vae$colinear
    ## the refinement ran AND took the held-intercept branch
    expect_gt(cd$nPhiTest, 0L)
    expect_gt(cd$nPhiClamp, 0L)
    ## and the bounds are actually respected in what comes back
    expect_lte(f$vae$zPop[2], log(2.0) + 1e-8)
    expect_lte(f$vae$zPop[3], log(25) + 1e-8)
  })

  test_that("a group larger than covSelectPhiMaxDim is skipped, not split", {
    skip_on_cran()
    ## Splitting an over-large correlated group would need an arbitrary
    ## tie-break on a near-tied graph, so the cap skips instead.  Without a test
    ## a regression could silently start evaluating the group, or splitting it,
    ## and nothing would fail.
    .three <- function() {
      ini({
        tka <- log(1.5); tcl <- log(2.7); tv <- log(31)
        eta.ka + eta.cl + eta.v ~ c(0.09,
                                    0.07, 0.09,
                                    0.07, 0.07, 0.09)
        add.sd <- 0.3
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    d <- .phiData(seed = 31L)
    ## a low join threshold so all three dims land in one component
    f <- suppressMessages(suppressWarnings(
      nlmixr2(.three, d, est = "vae",
              ## must run past klWarmup (50) or the refinement never gates on
              control = vaeControl(iters = 60L, itersBurnIn = 15L,
                                   calcTables = FALSE,
                                   covSelectPhiJoin = 0.2,
                                   covSelectPhiLeave = 0.1,
                                   covSelectPhiMaxDim = 2L))))
    cd <- f$vae$colinear
    ## the group formed...
    expect_gt(cd$nPhiPair, 0L)
    ## ...was refused for being too large...
    expect_gt(cd$nPhiSkipBig, 0L)
    ## ...and no move was scored, which is what "skipped, not split" means
    expect_identical(cd$nPhiTest, 0L)
    expect_identical(cd$nPhiMove, 0L)
  })

  test_that("the group-parallel refinement does not depend on the thread count", {
    skip_on_cran()
    ## The refinement runs one group per thread.  That is only sound because the
    ## correlated-dim components PARTITION the latent dims, so no two groups
    ## write the same cells, and because every group reads the residual snapshot
    ## taken before the region rather than live state another group is mutating.
    .old <- rxode2::rxCores()
    on.exit(rxode2::setRxThreads(.old), add = TRUE)
    d <- .phiData(seed = 31L)
    ctl <- vaeControl(iters = 60L, itersBurnIn = 15L, calcTables = FALSE,
                      covSelectPhiJoin = 0.4, covSelectPhiLeave = 0.3)
    rxode2::setRxThreads(1L)
    f1 <- suppressMessages(suppressWarnings(
      nlmixr2(.phiBlock, d, est = "vae", control = ctl)))
    rxode2::setRxThreads(8L)
    f8 <- suppressMessages(suppressWarnings(
      nlmixr2(.phiBlock, d, est = "vae", control = ctl)))
    f8b <- suppressMessages(suppressWarnings(
      nlmixr2(.phiBlock, d, est = "vae", control = ctl)))

    ## the mechanism must actually have done something, or this proves nothing
    expect_gt(f8$vae$colinear$nPhiMove, 0L)
    ## identical work, and identical decisions, whatever the thread count
    expect_identical(f1$vae$colinear$nPhiTest, f8$vae$colinear$nPhiTest)
    expect_identical(f1$vae$colinear$nPhiMove, f8$vae$colinear$nPhiMove)
    expect_identical(f1$vae$colinear$nPhiPair, f8$vae$colinear$nPhiPair)
    expect_identical(f1$vae$selected, f8$vae$selected)
    ## At a FIXED thread count the whole fit is bit-identical, which is the
    ## sharpest statement available: across DIFFERENT thread counts the VAE has
    ## pre-existing floating-point variation of order 1e-14 that predates this
    ## feature (it is present with covSelectColinear = FALSE), so the estimates
    ## are compared at a tolerance rather than for exact equality.
    expect_identical(f8$vae$beta, f8b$vae$beta)
    expect_identical(f8$vae$omega, f8b$vae$omega)
    expect_equal(f1$vae$beta, f8$vae$beta, tolerance = 1e-8)
  })

  test_that("hysteresis does not change a fit that has nothing to hold", {
    skip_on_cran()
    ## theo_sd offers one covariate group, so no cluster can bind and the whole
    ## colinear path is unreachable -- the selection must be bit-identical with
    ## the mechanism on and off
    skip_if_not_installed("nlmixr2data")
    d <- nlmixr2data::theo_sd
    mod <- .colinModel
    on <- suppressMessages(suppressWarnings(
      nlmixr2(mod, d, est = "vae", control = .ctl())))
    off <- suppressMessages(suppressWarnings(
      nlmixr2(mod, d, est = "vae", control = .ctl(covSelectColinear = FALSE))))
    expect_identical(on$vae$selected, off$vae$selected)
    expect_equal(on$vae$zPop, off$vae$zPop)
    expect_identical(on$vae$colinear$nColinTest, 0L)
  })
})
