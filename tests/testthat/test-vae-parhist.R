## The VAE parameter walk is central to the method: it is always captured (via
## the shared scale.h iteration-print machinery) and exposed as the standard
## $parHist. Training is RNG-based, so the seed (default 42) makes every fit
## reproducible.

.vaePhMod <- function() {
  ini({ lka <- log(1.8); lke <- log(0.086); lV <- log(32)
    eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03; add.err <- 0.7 })
  model({ ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
    d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
    cp <- central / V; cp ~ add(add.err) })
}

test_that("vaeControl seeds with 42 by default", {
  expect_equal(vaeControl()$seed, 42L)
})

test_that("est=vae captures the parameter walk and is reproducible", {
  skip_on_cran()
  ctl <- vaeControl(itersBurnIn = 5L, iters = 8L, klWarmup = 4L, gammaIter = 5L,
                    nGradStep = 2L, covariateSelection = FALSE, print = 0L)
  fit <- suppressMessages(suppressWarnings(
    nlmixr2(.vaePhMod(), nlmixr2data::theo_sd, est = "vae", control = ctl)))

  ## seed recorded on the fit
  expect_equal(fit$env$vae$seed, 42L)

  ## $parHist is the standard iter/objf + parameter data.frame with one row per
  ## captured iteration (burn-in + main)
  ph <- fit$parHist
  expect_s3_class(ph, "data.frame")
  expect_equal(nrow(ph), 5L + 8L)
  expect_true(all(c("iter", "objf", "lka", "lke", "lV", "add.err") %in% names(ph)))
  expect_false("type" %in% names(ph))            # .parHistCalc drops the type col
  expect_true(all(is.finite(ph$objf)))

  ## same seed -> identical walk
  fit2 <- suppressMessages(suppressWarnings(
    nlmixr2(.vaePhMod(), nlmixr2data::theo_sd, est = "vae", control = ctl)))
  expect_equal(fit$parHist, fit2$parHist)
})

test_that("est=vae does not perturb the caller's global RNG state", {
  skip_on_cran()
  ## the seed is set ONCE for the whole estimation (rxWithSeed in nlmixr2Est.vae)
  ## and the caller's global .Random.seed is restored on exit
  ctl <- vaeControl(itersBurnIn = 3L, iters = 4L, klWarmup = 2L, gammaIter = 3L,
                    nGradStep = 2L, covariateSelection = FALSE, print = 0L)
  set.seed(123); want <- runif(3)
  set.seed(123)
  invisible(suppressMessages(suppressWarnings(
    nlmixr2(.vaePhMod(), nlmixr2data::theo_sd, est = "vae", control = ctl))))
  got <- runif(3)
  expect_equal(got, want)
})
