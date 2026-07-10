## VAE C++ inner-likelihood driver: set up the FOCEi inner problem once
## (.vaeInnerSetup -> vaeInnerSetup_) then evaluate the per-subject objective and
## eta-gradient through likInner0/lpInner in parallel (.vaeInnerEval ->
## vaeInnerLik). Must match the R-exported likInner/foceiInnerLp, and for a
## mixture model must give distinct per-component objectives (hard assignment).

test_that("vae inner driver matches likInner/foceiInnerLp (theophylline)", {
  theo <- function() {
    ini({ lka <- log(1.8); lke <- log(0.086); lV <- log(32)
      eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03; add.err <- 0.7 })
    model({ ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
      d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
      cp <- central / V; cp ~ add(add.err) })
  }
  ui <- rxode2::assertRxUi(theo)
  ctl <- vaeControl()
  N <- length(unique(nlmixr2data::theo_sd$ID))
  set.seed(1); etaMat <- matrix(rnorm(N * 3, 0, 0.1), N, 3)

  .vaeInnerSetup(ui, nlmixr2data::theo_sd, etaMat, ctl)
  on.exit(.vaeInnerFree(), add = TRUE)
  r <- .vaeInnerEval(etaMat, ctl, grad = TRUE)
  expect_equal(length(r$obj), N)
  ## match the single-subject R-exported inner functions
  expect_lt(abs(r$obj[1] - likInner(etaMat[1, ], 1L)), 1e-6)
  expect_lt(max(abs(r$lp[1, ] - foceiInnerLp(etaMat[1, ], 1L))), 1e-6)
  expect_lt(abs(r$obj[5] - likInner(etaMat[5, ], 5L)), 1e-6)
})

test_that("vae inner driver selects mixture components per id", {
  mixmod <- function() {
    ini({ lka <- log(1.8); lke1 <- log(0.15); lke2 <- log(0.04); lV <- log(32); p1 <- 0.6
      eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03; add.err <- 0.7 })
    model({ ka <- exp(lka + eta.ka)
      ke <- mix(exp(lke1 + eta.ke), p1, exp(lke2 + eta.ke)); V <- exp(lV + eta.V)
      d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
      cp <- central / V; cp ~ add(add.err) })
  }
  ui <- rxode2::assertRxUi(mixmod)
  expect_equal(ui$saemNMix, 2L)
  ctl <- vaeControl()
  N <- length(unique(nlmixr2data::theo_sd$ID)); nMix <- ui$saemNMix
  etaSetup <- matrix(0, N, 3)

  .vaeInnerSetup(ui, nlmixr2data::theo_sd, etaSetup, ctl)
  on.exit(.vaeInnerFree(), add = TRUE)
  ## nSub*nMix ids, component-major (encoder etas repeated per component)
  etaEval <- do.call(rbind, rep(list(etaSetup), nMix))
  r <- .vaeInnerEval(etaEval, ctl)
  expect_equal(length(r$obj), N * nMix)
  o1 <- r$obj[1:N]; o2 <- r$obj[(N + 1):(2 * N)]
  ## the two components give distinct per-subject objectives (mixture selection)
  expect_true(all(abs(o1 - o2) > 1e-6))
  ## hard-assignment mixnum is a valid component per subject
  mixnum <- ifelse(o1 <= o2, 1L, 2L)
  expect_true(all(mixnum %in% c(1L, 2L)))
})
