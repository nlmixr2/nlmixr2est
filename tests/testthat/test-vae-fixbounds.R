## Fixed parameters, fixed omega, non-mu-referenced (free) etas, and parameter
## bounds for est="vae". Structural/error thetas are estimated on their natural
## scale; the M-step holds fixed omega, forces free-eta centers to 0, and clamps
## population estimates to their [lower, upper] bounds (constrained estimate).

.vaeFbCtl <- function()
  vaeControl(itersBurnIn = 20L, iters = 50L, klWarmup = 15L, gammaIter = 35L,
             nGradStep = 3L, covariateSelection = FALSE, returnVae = TRUE, seed = 1L)

test_that("vae holds a fixed omega and a fixed mu-referenced theta", {
  skip_on_cran()
  mod <- function() {
    ini({ lka <- log(1.8); lke <- log(0.086); lV <- fixed(log(32))
      eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ fixed(0.05); add.err <- 0.7 })
    model({ ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
      d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
      cp <- central / V; cp ~ add(add.err) })
  }
  fit <- suppressMessages(suppressWarnings(
    nlmixr2(mod, nlmixr2data::theo_sd, est = "vae", control = .vaeFbCtl())))
  ## fixed omega held exactly at its initial value
  expect_equal(unname(fit$omega[3]), 0.05)
  ## the fixed mu-referenced theta (lV) drops to a free eta centered at 0
  expect_equal(unname(fit$zPop[3]), 0)
  ## the free (non-mu-ref) parameters are still estimated
  expect_true(is.finite(fit$zPop[1]) && is.finite(fit$omega[1]) && fit$omega[1] > 0)
})

test_that("vae clamps a structural parameter to its bound", {
  skip_on_cran()
  ## true lka ~ log(1.8)=0.588; a lower bound above it forces the estimate there
  mod <- function() {
    ini({ lka <- c(log(2.5), log(3), log(8)); lke <- log(0.086); lV <- log(32)
      eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03; add.err <- 0.7 })
    model({ ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
      d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
      cp <- central / V; cp ~ add(add.err) })
  }
  fit <- suppressMessages(suppressWarnings(
    nlmixr2(mod, nlmixr2data::theo_sd, est = "vae", control = .vaeFbCtl())))
  expect_gte(fit$zPop[1], log(2.5) - 1e-6)
  expect_lte(fit$zPop[1], log(8) + 1e-6)
  expect_lt(abs(fit$zPop[1] - log(2.5)), 1e-3)   # sits at the active lower bound
})

test_that("vae clamps an error parameter to its bound", {
  skip_on_cran()
  set.seed(3)
  sm <- function() {
    ini({ lke <- log(0.1); lV <- log(32) })
    model({ ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
      d/dt(central) = -ke * central; cp <- central / V })
  }
  ev <- rxode2::et(amt = 320, cmt = "central") %>%
    rxode2::et(seq(1, 24, length.out = 8)) %>% rxode2::et(id = 1:30)
  bs <- rxode2::rxSolve(sm, ev, params = c(lke = log(0.1), lV = log(32)),
                        omega = lotri::lotri(eta.ke ~ 0.02, eta.V ~ 0.02))
  bs <- as.data.frame(bs)[, c("id", "time", "cp")]; names(bs) <- c("ID", "TIME", "DV")
  bs$DV <- bs$DV + stats::rnorm(nrow(bs), 0, 0.15 * abs(bs$DV))
  dose <- data.frame(ID = unique(bs$ID), TIME = 0, DV = 0, EVID = 1, AMT = 320, CMT = 1)
  bs$EVID <- 0; bs$AMT <- 0; bs$CMT <- 1
  bs <- rbind(dose, bs); bs <- bs[order(bs$ID, bs$TIME, -bs$EVID), ]
  pm <- function() {
    ini({ lke <- log(0.1); lV <- log(32); eta.ke ~ 0.02; eta.V ~ 0.02
      prop.err <- c(0, 0.02, 0.05) })
    model({ ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
      d/dt(central) = -ke * central; cp <- central / V; cp ~ prop(prop.err) })
  }
  fit <- suppressMessages(suppressWarnings(
    nlmixr2(pm, bs, est = "vae", control = .vaeFbCtl())))
  expect_lte(fit$a[["prop.err"]], 0.05 + 1e-6)          # clamped to the upper bound
  expect_gt(fit$a[["prop.err"]], 0.04)
})
