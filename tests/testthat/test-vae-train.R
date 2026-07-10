## Phase 4 (Milestone A): end-to-end VAE training on theophylline WITHOUT
## covariate selection. A short schedule (kept small for CI) must run the full
## burn-in -> KL-anneal -> smoothing loop and land the fixed effects near the
## paper's Table 1 VAE column. The full-schedule parity check lives in the
## vignette; this is a fast smoke + sanity-of-direction test.

test_that("vae trains end to end and approaches theophylline Table 1", {
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
  ctl <- vaeControl(itersBurnIn = 30L, klWarmup = 30L, gammaIter = 60L,
                    iters = 80L, hiddenDim = 25L, seed = 1L)
  prep <- .vaeDataPrep(ui, nlmixr2data::theo_sd)
  ## training drives the FOCEi inner likelihood directly (set up once)
  innerEnv <- .vaeInnerSetup(ui, nlmixr2data::theo_sd, matrix(0, prep$N, prep$zDim), ctl)
  on.exit(.vaeInnerFree(), add = TRUE)
  fit <- .vaeTrain(prep, innerEnv, ctl)

  ka <- exp(fit$zPop[1]); ke <- exp(fit$zPop[2]); V <- exp(fit$zPop[3])
  ## ELBO trace should have descended and estimates be finite
  expect_true(all(is.finite(fit$zPop)))
  expect_true(all(is.finite(fit$omega)) && all(fit$omega > 0))
  expect_true(is.finite(fit$a) && fit$a > 0)
  expect_lt(fit$elboTrace[length(fit$elboTrace)], fit$elboTrace[1])

  ## fixed effects within a loose band of the paper VAE column (short schedule)
  expect_lt(abs(ka - 1.63) / 1.63, 0.10)
  expect_lt(abs(ke - 0.0867) / 0.0867, 0.10)
  expect_lt(abs(V - 31.97) / 31.97, 0.10)
  expect_lt(abs(fit$a - 0.71) / 0.71, 0.10)
})
