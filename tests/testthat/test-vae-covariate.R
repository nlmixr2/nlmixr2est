## Phase 5 (Milestone B): BICc-ELBO covariate selection on theophylline. The VAE
## must auto-discover WT and select it on ka and V (not ke), matching the paper,
## and the WT effect on ka must pull omega_ka down from its no-covariate value.
## (Slow: runs a moderate training schedule; skipped on CRAN.)

nmTest({
  test_that("vae selects WT on ka and V (not ke) for theophylline", {
    skip_on_cran()
    theo <- function() {
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
    ui <- rxode2::assertRxUi(theo)
    prep <- .vaeDataPrep(ui, nlmixr2data::theo_sd)
    expect_equal(prep$covNames, "WT")
    expect_equal(prep$covType, "continuous")

    ctl <- vaeControl(itersBurnIn = 80L, klWarmup = 40L, gammaIter = 120L,
                      iters = 160L, hiddenDim = 25L, seed = 1L, covariateSelection = TRUE,
                      print = 0L)
    innerEnv <- .vaeInnerSetup(ui, nlmixr2data::theo_sd,
                               matrix(0, prep$N, prep$zDim), ctl)
    on.exit(.vaeInnerFree(), add = TRUE)
    fit <- rxode2::rxWithSeed(1L, .vaeTrain(prep, innerEnv, ctl))

    ## rows = params (ka, ke, V), single column WT
    expect_true(fit$selected[1, 1])   # WT -> ka
    expect_false(fit$selected[2, 1])  # WT -> ke NOT selected
    expect_true(fit$selected[3, 1])   # WT -> V
    ## WT->ka effect is large and positive (paper ~2.55)
    expect_gt(fit$beta[1, 1], 1.5)
    ## omega_ka pulled down toward the with-covariate value (< no-covariate ~0.61)
    expect_lt(sqrt(fit$omega[1]), 0.58)
    ## fixed effects sane
    expect_lt(abs(exp(fit$zPop[2]) - 0.0867) / 0.0867, 0.05)   # ke
    expect_lt(abs(exp(fit$zPop[3]) - 31.97) / 31.97, 0.10)     # V
  })

  ## End-to-end regression for the two covariate-output bugs, exercised together
  ## via a fixed residual parameter (literalFix=TRUE default), which is what
  ## triggered the parFixedDf drop:
  ##   Bug 1 -- the selected covariate coefficients (beta_*) must appear in the
  ##            population-parameter table, not just in $theta/$cov.
  ##   Bug 2 -- covariate-bearing mu-parameters must back-transform (exp) rather
  ##            than print the raw log-scale estimate.
  test_that("est=vae keeps covariate betas and back-transforms with a fixed parameter", {
    skip_on_cran()
    theoFix <- function() {
      ini({
        lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03
        add.err <- fix(0.7)
      })
      model({
        ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot
        d/dt(central) = ka * depot - ke * central
        cp <- central / V
        cp ~ add(add.err)
      })
    }
    ctl <- vaeControl(itersBurnIn = 80L, klWarmup = 40L, gammaIter = 120L,
                      iters = 160L, hiddenDim = 25L, seed = 1L,
                      covariateSelection = TRUE, print = 0L)
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(theoFix, nlmixr2data::theo_sd, est = "vae", control = ctl)))

    pf <- fit$parFixedDf
    ## Bug 1: covariate coefficients present in the population-parameter table
    ## (WT selected on ka and V for theophylline)
    expect_true(all(c("beta_lka_WT", "beta_lV_WT") %in% rownames(pf)))
    expect_true(all(rownames(pf) %in% rownames(fit$cov) |
                      rownames(pf) == "add.err"))

    ## Bug 2: the covariate-bearing mu-parameters back-transform (exp), so the
    ## Back-transformed value differs from the raw log-scale Estimate
    expect_equal(pf["lka", "Back-transformed"], exp(pf["lka", "Estimate"]),
                 tolerance = 1e-6)
    expect_equal(pf["lV", "Back-transformed"], exp(pf["lV", "Estimate"]),
                 tolerance = 1e-6)
    ## the covariate-free lke also back-transforms
    expect_equal(pf["lke", "Back-transformed"], exp(pf["lke", "Estimate"]),
                 tolerance = 1e-6)
    ## the fixed residual parameter is on the natural scale (no exp)
    expect_equal(pf["add.err", "Estimate"], 0.7, tolerance = 1e-6)
    ## covariate coefficients are reported raw (not back-transformed)
    expect_equal(pf["beta_lka_WT", "Back-transformed"],
                 pf["beta_lka_WT", "Estimate"], tolerance = 1e-6)
  })

  ## The L0-penalty warmup ramp (covSelectAlpha) is a distinct step in the
  ## iteration table -- it must be labeled "CovSel ramp" on the ramp iterations,
  ## which are exactly the ones that evaluate and print the ELBO objective.
  test_that("est=vae shows the CovSel ramp step in the iteration print", {
    skip_on_cran()
    theo <- function() {
      ini({ lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03; add.err <- 0.7 })
      model({ ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ctl <- vaeControl(itersBurnIn = 2L, iters = 6L, klWarmup = 4L, gammaIter = 5L,
                      nGradStep = 2L, covariateSelection = TRUE, covSelectAlpha = 2,
                      print = 1L)
    tc <- textConnection("msgs", "w", local = TRUE)
    sink(tc, type = "message")
    out <- capture.output(
      suppressWarnings(nlmixr2(theo, nlmixr2data::theo_sd, est = "vae", control = ctl)),
      type = "output")
    sink(type = "message"); close(tc)
    allout <- c(out, msgs)
    ## the ramp iterations (it < klWarmup) are labeled, and the legend documents it
    expect_true(any(grepl("CovSel ramp", allout)))
    expect_true(any(grepl("covariate-selection L0-penalty warmup", allout)))
  })
})
