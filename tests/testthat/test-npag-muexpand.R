## saem-style mu-expansion for npag: a non-mu-referenced structural fixed-effect
## theta (err==NA, no eta -- e.g. `ke <- exp(tke)`) is not a grid dimension and
## feeds the ODE, so npag cannot estimate it directly.  The mu-expansion injects a
## pseudo-eta (`ke <- exp(tke + eta.tke)`) so the parameter becomes a grid
## dimension estimated by the support-point distribution (same convention as every
## npag theta: the reference theta stays put and the location lives in the
## support).  Real fit -> weekly slow batch.

nmTest({
  # deliberately-wrong ke start (0.30); theo's ke is ~0.08-0.10
  .mod <- function() {
    ini({ tka <- log(1.5); tv <- log(32); tke <- log(0.30)
      eta.ka ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); ke <- exp(tke)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - ke * center
      cp <- center / v; cp ~ add(add.sd) })
  }
  .ctl <- function(muExpand = TRUE)
    npagControl(points = 128L, cycles = 30L, gammaOptimize = FALSE, seed = 1L,
                calcTables = FALSE, muExpand = muExpand)

  test_that("est='npag' mu-expansion recovers a non-mu structural theta as a fixed effect", {
    f <- nlmixr2(.mod, nlmixr2data::theo_sd, est = "npag", control = .ctl(TRUE))
    expect_true("eta.tke" %in% rownames(f$omega))
    # the injected eta's support-mean is folded into the theta and the random effect
    # is collapsed, so the REPORTED theta is the estimate (recovers theo ke from the
    # deliberately-wrong 0.30 start) and there is no BSV on the fixed effect.
    expect_true(exp(as.numeric(f$theta[["tke"]])) > 0.05 &&
                exp(as.numeric(f$theta[["tke"]])) < 0.15)
    expect_lt(f$omega["eta.tke", "eta.tke"], 1e-4)     # collapsed: fixed effect, no BSV
    .sp <- f$env$npagSupport
    expect_true(all(abs(.sp[, ncol(.sp)]) < 1e-8))     # injected eta support zeroed
  })

  test_that("est='npag' muExpand=FALSE leaves the non-mu structural theta held", {
    f <- nlmixr2(.mod, nlmixr2data::theo_sd, est = "npag", control = .ctl(FALSE))
    expect_false("eta.tke" %in% rownames(f$omega))
    expect_equal(as.numeric(f$theta[["tke"]]), log(0.30))   # held at ini
  })

  # fixture 2: a non-mu-referenced eta (enters non-additively, so rxode2 does not
  # pair it with a theta) needs NO mu-expansion -- the npag box covers every eta, so
  # it is already a grid dimension estimated as a pure random effect (population mean
  # ~0).  Confirm it fits and its support distribution is recovered.
  .modEta <- function() {
    ini({ tka <- log(1.5); tv <- log(32); tke <- log(0.1)
      eta.ka ~ 0.3; eta.v ~ 0.1; eta.ke ~ 0.1; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); ke <- exp(tke) * (1 + eta.ke)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - ke * center
      cp <- center / v; cp ~ add(add.sd) })
  }

  test_that("est='npag' natively estimates a non-mu-referenced eta as a grid dimension", {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.modEta))
    expect_false("eta.ke" %in% .ui$muRefDataFrame$eta)   # genuinely non-mu
    f <- nlmixr2(.modEta, nlmixr2data::theo_sd, est = "npag",
                 control = npagControl(points = 64L, cycles = 10L, gammaOptimize = FALSE,
                                       seed = 1L, calcTables = FALSE))
    expect_true("eta.ke" %in% rownames(f$omega))          # it is a support dimension
    .sp <- f$env$npagSupport; .wt <- f$env$npagWeights
    .j <- ncol(.sp)                                        # eta.ke is the last eta
    expect_true(abs(sum(.wt * .sp[, .j])) < 0.5)           # pure RE: support mean ~ 0
    expect_gt(f$omega["eta.ke", "eta.ke"], 0)              # a BSV was estimated
  })
})
