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
})
