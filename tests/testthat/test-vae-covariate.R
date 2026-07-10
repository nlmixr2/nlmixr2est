## Phase 5 (Milestone B): BICc-ELBO covariate selection on theophylline. The VAE
## must auto-discover WT and select it on ka and V (not ke), matching the paper,
## and the WT effect on ka must pull omega_ka down from its no-covariate value.
## (Slow: runs a moderate training schedule; skipped on CRAN.)

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
  am <- .vaeDecoderModel(ui)

  ctl <- vaeControl(itersBurnIn = 80L, klWarmup = 40L, gammaIter = 120L,
                    iters = 160L, hiddenDim = 25L, seed = 1L, covariateSelection = TRUE)
  fit <- .vaeTrain(prep, am, ctl)

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
