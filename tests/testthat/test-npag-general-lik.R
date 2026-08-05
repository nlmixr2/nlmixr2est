## Generalized (non-normal) likelihood support for NPAG.  The npag Psi sums the
## inner per-observation llikObs, which for a non-normal endpoint is exactly the
## user's log-likelihood -- so the nonparametric objective is already correct.
## The residual/likelihood parameters (iniDf$err non-NA, e.g. a t-distribution's
## df) are estimated with the same frozen-ODE bounded step as the residual params;
## the ODE freeze is valid because those parameters feed only the post-solve f/r,
## never the state derivatives.  npb (Gibbs) still rejects non-normal endpoints.
## Real fit -> weekly slow batch.

nmTest({
  .tMod <- function() {
    ini({ tka <- log(1.5); tv <- log(32); tke <- fix(log(0.08))   # fixed: focus on the likelihood
      eta.ka ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7; nu <- 8 })
    model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); ke <- exp(tke)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - ke * center
      cp <- center / v
      cp ~ add(add.sd) + dt(nu) })          # Student-t residual: a general likelihood
  }

  test_that("est='npag' fits a generalized (Student-t) likelihood and estimates its params", {
    f <- nlmixr2(.tMod, nlmixr2data::theo_sd, est = "npag",
                 control = npagControl(points = 48L, cycles = 8L, seed = 1L))
    expect_s3_class(f, "nlmixr2FitData")
    expect_true(is.finite(as.numeric(f$objf)))
    # the t-distribution parameters (err-tagged: add / t) are estimated, not held
    expect_true(as.numeric(f$theta[["add.sd"]]) > 0)
    expect_true(as.numeric(f$theta[["nu"]]) > 0)
    expect_false(isTRUE(all.equal(as.numeric(f$theta[["nu"]]), 8)))   # moved off the start
    # gamma is meaningless for a non-normal endpoint (r == 1) -> forced off
    expect_false(isTRUE(f$control$npGammaOptimize))
    # only err-tagged params are optimized (tke is fixed), so the ODE freeze is valid
    expect_true(isTRUE(f$control$npResidFreeze))
  })

  test_that("est='npb' also fits a generalized (Student-t) likelihood", {
    # npb sums the same per-observation llikObs as npag, so the Gibbs sweep handles
    # a non-normal endpoint too (it just samples the mixing distribution rather than
    # optimizing an err-tagged likelihood parameter).
    f <- nlmixr2(.tMod, nlmixr2data::theo_sd, est = "npb",
                 control = npbControl(points = 20L, burnin = 15L, nsamp = 15L, seed = 1L))
    expect_s3_class(f, "nlmixr2FitData")
    expect_true(is.finite(as.numeric(f$objf)))
    expect_equal(ncol(f$env$npbSupport), 2L)   # eta.ka, eta.v support dimensions
  })

  ## ---- MULTIPLE endpoints (#838) ------------------------------------------------
  ## Everything above is SINGLE-endpoint, which is exactly the case #838 could never
  ## reach: with one endpoint rx_yj_ is a constant, so reading the endpoint's
  ## distribution before calc_lhs had evaluated the row still gave the right answer.
  ##
  ## npag/npb have no defect of their own here -- npEvalCondLik just sums likInner0's
  ## per-observation llikObs -- but that is precisely why they inherited #838: one row
  ## per subject entered Psi scored as NORMAL, i.e. with its log-density treated as a
  ## prediction of DV against a variance forced to 1, contributing -0.5*(logDens - DV)^2.
  ## Measured on the model below at points=32/cycles=4: objf 162145.67 before the fix
  ## against -5326.82 after, while the all-Gaussian twin was bit-identical across both
  ## builds (133.907076101) -- so the swing is confined to the general-likelihood path.
  .pkpdLL <- function() {
    ini({ tka <- 0.5; tcl <- -3.2; tv <- -0.7; tec50 <- 2; tkout <- -2; te0 <- 4.6
      eta.cl ~ 0.09; add.sd <- 0.4; pdadd.sd <- 2 })
    model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      ec50 <- exp(tec50); kout <- exp(tkout); e0 <- exp(te0)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      effect(0) <- e0
      d/dt(effect) <- kout * (e0 * (1 - cp / (ec50 + cp)) - effect)
      cp ~ add(add.sd) | center
      ll(effect) ~ -0.5 * log(2 * pi) - log(pdadd.sd) -
        0.5 * ((DV - effect) / pdadd.sd)^2 })
  }
  .pkpdGauss <- function() {
    ini({ tka <- 0.5; tcl <- -3.2; tv <- -0.7; tec50 <- 2; tkout <- -2; te0 <- 4.6
      eta.cl ~ 0.09; add.sd <- 0.4; pdadd.sd <- 2 })
    model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      ec50 <- exp(tec50); kout <- exp(tkout); e0 <- exp(te0)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      effect(0) <- e0
      d/dt(effect) <- kout * (e0 * (1 - cp / (ec50 + cp)) - effect)
      cp ~ add(add.sd) | center
      effect ~ add(pdadd.sd) | effect })
  }
  .npPkpdData <- function() {
    rxode2::rxSetSeed(1)
    .ev <- rxode2::et(amt = 100, cmt = "depot", id = 1:12)
    .ev <- rxode2::et(.ev, seq(0.5, 24, by = 3), cmt = "center")
    .ev <- rxode2::et(.ev, seq(0.5, 24, by = 3), cmt = "effect")
    .d <- as.data.frame(rxode2::rxSolve(.pkpdGauss(), .ev, addDosing = TRUE))
    .dose <- .d[.d$evid != 0, c("id", "time", "CMT", "amt", "evid")]
    .dose$dv <- NA_real_
    names(.dose)[names(.dose) == "CMT"] <- "cmt"
    .obs <- .d[.d$evid == 0, c("id", "time", "CMT", "sim")]
    .obs$amt <- 0; .obs$evid <- 0
    names(.obs)[names(.obs) == "CMT"] <- "cmt"
    names(.obs)[names(.obs) == "sim"] <- "dv"
    .dat <- rbind(.dose[, c("id", "time", "dv", "cmt", "amt", "evid")],
                  .obs[, c("id", "time", "dv", "cmt", "amt", "evid")])
    .dat[order(.dat$id, .dat$time, -.dat$evid), ]
  }
  ## A Gaussian (or dnorm()) twin does NOT work as a self-validating reference here,
  ## unlike the focei test in test-focei-ll-fast-grad-fit.R.  Both were tried: npag's
  ## objective is a NONPARAMETRIC mixture likelihood, so the two models converge to
  ## different support sets and the gap is not a constant.  Measured at
  ## points=24/cycles=3: hand-written ll() -5418.12, add() twin 90.70, the same density
  ## written as `add(pdadd.sd) + dnorm()` 925.30.  Do not "fix" this into a twin
  ## comparison -- assert the failure mode instead, which is orders of magnitude away.
  .npObjfOk <- function(o) {
    expect_true(is.finite(o))
    # the defect makes this a large POSITIVE data-sized excess (+1.6e5 measured);
    # a collapsed/degenerate fit would run off the other way
    expect_lt(o, -1000)
    expect_gt(o, -5e4)
  }

  test_that("est='npag' scores a general likelihood alongside a second endpoint (#838)", {
    skip_on_cran()
    .dat <- .npPkpdData()
    # the mechanism under test: two endpoints, at least one of them non-normal
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.pkpdLL))
    expect_equal(nrow(.ui$predDf), 2L)
    expect_false(all(as.character(.ui$predDfFocei$distribution) %in% c("norm", "dnorm")))
    f <- suppressWarnings(nlmixr2(.pkpdLL, .dat, est = "npag",
                                  control = npagControl(points = 24L, cycles = 3L,
                                                        seed = 1L, print = 0L)))
    expect_s3_class(f, "nlmixr2FitData")
    .npObjfOk(as.numeric(f$objf))
  })

  test_that("est='npb' scores a general likelihood alongside a second endpoint (#838)", {
    skip_on_cran()
    .dat <- .npPkpdData()
    f <- suppressWarnings(nlmixr2(.pkpdLL, .dat, est = "npb",
                                  control = npbControl(points = 16L, burnin = 8L,
                                                       nsamp = 8L, seed = 1L, print = 0L)))
    expect_s3_class(f, "nlmixr2FitData")
    .npObjfOk(as.numeric(f$objf))
  })
})
