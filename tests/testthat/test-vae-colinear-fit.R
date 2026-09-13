## The colinearity cluster reaching the C++ covariate M-step: the hysteresis
## pass (which keeps last iteration's incumbent unless a cluster mate beats it
## by a covariate's L0 cost) and the near-tie record.  Both live behind ONE gate
## -- R ships the cluster vector only when .vaeClusterBinds() is TRUE -- so the
## pair of tests here drives the gate open and shut.  Fit-based: .slowBatches.

nmTest({

  ## theo_sd plus a near-duplicate of WT.  Both are subject-constant, so both
  ## survive the time-varying screen and reach the search; they correlate well
  ## past the default 0.9 cut, so they land in one cluster.
  .theoColinear <- function(seed = 7L) {
    d <- nlmixr2data::theo_sd
    .id <- unique(d$ID)
    .wt <- vapply(.id, function(i) d$WT[d$ID == i][1], numeric(1))
    .lbm <- rxode2::rxWithSeed(seed,
                               .wt * 0.8 + stats::rnorm(length(.wt), sd = 0.4))
    d$LBM <- .lbm[match(d$ID, .id)]
    list(data = d, cor = abs(stats::cor(.wt, .lbm)))
  }

  .theo <- function() {
    ini({
      lka <- log(1.8); lke <- log(0.086); lV <- log(32)
      eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03
      add.err <- 0.7
    })
    model({
      ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
      d/dt(depot) = -ka * depot
      d/dt(central) = ka * depot - ke * central
      cp <- central / V
      cp ~ add(add.err)
    })
  }

  .runVae <- function(d, cut) {
    ctl <- vaeControl(itersBurnIn = 40L, klWarmup = 20L, gammaIter = 60L,
                      iters = 90L, hiddenDim = 15L, seed = 1L,
                      covariateSelection = TRUE, print = 0L,
                      shapes = "power", covCenterType = "mean",
                      covSelectColinearCut = cut)
    ui <- rxode2::assertRxUi(.theo)
    prep <- .vaeDataPrep(ui, d, ctl)
    innerEnv <- .vaeInnerSetup(ui, d, matrix(0, prep$N, prep$zDim), ctl)
    on.exit(.vaeInnerFree(), add = TRUE)
    list(prep = prep, fit = rxode2::rxWithSeed(1L, .vaeTrain(prep, innerEnv, ctl)))
  }

  test_that("a colinear pair opens the gate: hysteresis runs, near ties report", {
    skip_on_cran()
    .d <- .theoColinear()
    ## pin the premise the default cut keys on, so a change in the fixture is
    ## not mistaken for a change in the mechanism
    expect_gt(.d$cor, 0.9)
    .r <- .runVae(.d$data, cut = .vaeColinearCut)
    expect_true(all(c("WT_power", "LBM_power") %in% .r$prep$covNames))
    ## The hysteresis counter is reported whenever the gate is open.  It counts
    ## INTERVENTIONS, so it stays 0 unless the search actually flips between two
    ## mates -- measured not to happen here even at cor = 1, because the leaf
    ## tie-break is deterministic.  The decision itself is tested directly, in
    ## test-vae-colinear.R, through vaeClusterSwapOnly_().
    expect_true(is.integer(.r$fit$nCovHysteresis))
    expect_gte(.r$fit$nCovHysteresis, 0L)
    ## near ties: a mate is never the column itself, both name real covariates,
    ## and the mate is by construction no better than what was chosen
    nt <- .r$fit$covNearTie
    expect_gt(nrow(nt), 0L)
    expect_true(all(nt$covariate != nt$mate))
    expect_true(all(c(nt$covariate, nt$mate) %in% .r$prep$covNames))
    expect_true(all(nt$delta >= 0))
    expect_true(all(nt$eta %in% .r$prep$etaNames))
  })

  test_that("covSelectColinearCut = 1 shuts the gate", {
    skip_on_cran()
    .d <- .theoColinear()
    ## only an exact duplicate clusters at 1, so .vaeClusterBinds() is FALSE,
    ## R ships no cluster vector and neither mechanism can run
    .r <- .runVae(.d$data, cut = 1)
    expect_identical(.r$fit$nCovHysteresis, 0L)
    expect_identical(nrow(.r$fit$covNearTie), 0L)
  })

  ## -----------------------------------------------------------------------
  ## Cross-parameter refinement: the second, latent-dim grouping pass.  It is
  ## independent of the covariate cluster above -- it groups latent DIMS by the
  ## empirical correlation of the posterior means, and it only acts under a
  ## correlated omega.  These need real fits (the pass lives inside training).
  ## -----------------------------------------------------------------------

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
    ## the gate opened, because the model declares the correlation
    expect_true(f$vae$omOff)
    ## groups formed and moves were evaluated -- "the mechanism ran"
    expect_gt(f$vae$nPhiPair, 0L)
    expect_gt(f$vae$nPhiTest, 0L)
    ## no advisory: the refinement was not skipped, so there is nothing to advise
    expect_false(any(grepl("declare an omega block", f$runInfo)))
  })

  test_that("a diagonal omega detects the groups but skips the refinement", {
    skip_on_cran()
    f <- suppressMessages(suppressWarnings(
      nlmixr2(.phiDiag, .phiData(), est = "vae", control = .phiCtl())))
    expect_false(f$vae$omOff)
    ## the dims ARE correlated and that is reported...
    expect_gt(f$vae$nPhiPair, 0L)
    ## ...but with a diagonal omega the objective is separable across dims, so
    ## not one move is scored -- this is a proof of the gate, not of an empty
    ## result
    expect_identical(f$vae$nPhiTest, 0L)
    expect_identical(f$vae$nPhiMove, 0L)
    ## and the modeler is told what would enable it
    expect_match(f$runInfo, "declare an omega block", all = FALSE)
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
      expect_gt(f$vae$nPhiPair, 0L, label = src)
      expect_gt(f$vae$nPhiTest, 0L, label = src)
    }
  })

  test_that("a wider sticky band never joins fewer pairs", {
    skip_on_cran()
    ## covSelectPhiLeave below covSelectPhiJoin is what keeps a pair joined when
    ## its correlation dips into the gap.  Widening the band can therefore only
    ## ever ADD joined pair-iterations, never remove them, whatever the
    ## correlation trajectory happens to be.  $vae$phiPairOn exposes the
    ## adjacency for anyone constructing a sharper fixture later.
    d <- .phiData(seed = 31L)
    .pairs <- function(lv) {
      f <- suppressMessages(suppressWarnings(
        nlmixr2(.phiBlock, d, est = "vae",
                control = vaeControl(iters = 80L, itersBurnIn = 15L,
                                     calcTables = FALSE,
                                     covSelectPhiCor = "mu",
                                     covSelectPhiJoin = 0.78,
                                     covSelectPhiLeave = lv))))
      list(n = f$vae$nPhiPair, on = f$vae$phiPairOn)
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
    ## the refinement ran AND took the held-intercept branch
    expect_gt(f$vae$nPhiTest, 0L)
    expect_gt(f$vae$nPhiClamp, 0L)
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
    ## the group formed...
    expect_gt(f$vae$nPhiPair, 0L)
    ## ...was refused for being too large...
    expect_gt(f$vae$nPhiSkipBig, 0L)
    ## ...and no move was scored, which is what "skipped, not split" means
    expect_identical(f$vae$nPhiTest, 0L)
    expect_identical(f$vae$nPhiMove, 0L)
  })
})
