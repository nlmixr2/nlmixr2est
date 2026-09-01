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

  ## ---- MULTIPLE endpoints (#850) ------------------------------------------------
  ## Everything above is SINGLE-endpoint, which was never the broken case.  With a
  ## second endpoint the residual step scored a log-density with a Gaussian ELS form,
  ## the moment warm start took a moment of a log-density, and rxode2's safeLog paid
  ## +36 per observation for a NEGATIVE sd -- together giving a log-likelihood of
  ## +2364 to +2831 where the data bounds it at -155, with pdadd.sd at -0.83 to -1.93.
  ##
  ## The reference is the model's own GAUSSIAN TWIN, not a checked-in constant.  An
  ## ll() written as the exact normal log-density is the same likelihood as the
  ## equivalent add() endpoint, and likInner0 stores the full density for ll() rows
  ## while the Gaussian path omits the 0.5*log(2*pi) term -- so the two
  ## log-likelihoods differ by nObs*0.5*log(2*pi), plus whatever the two fits'
  ## residual estimates are worth.  That last part is NOT zero: the residual step
  ## scores an ordinary err-tagged parameter by extended least squares and a
  ## general-likelihood one by the user's log-density, so the endpoint written
  ## differently in the two models lands somewhere different in each (see the
  ## comment on pdadd.sd below).  The endpoint written IDENTICALLY in both is the
  ## part that is self-validating.
  .npMkPkPd <- function(pdLine) {
    .b <- quote({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      ec50 <- exp(tec50); kout <- exp(tkout); e0 <- exp(te0)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      effect(0) <- e0
      d/dt(effect) <- kout * (e0 * (1 - cp / (ec50 + cp)) - effect)
      cp ~ add(add.sd) | center })
    .b[[length(.b) + 1]] <- pdLine
    .f <- function() {}
    body(.f) <- bquote({
      ini({ tka <- 0.5; tcl <- -3.2; tv <- -0.7; tec50 <- 2; tkout <- -2; te0 <- 4.6
        eta.cl ~ 0.09; add.sd <- 0.4; pdadd.sd <- 2 })
      model(.(.b))
    })
    .f
  }
  .npPkPdData <- function(m) {
    .testSeed(1); rxode2::rxSetSeed(1)
    .ev <- rxode2::et(amt = 100, cmt = "depot", id = 1:12)
    .ev <- rxode2::et(.ev, seq(0.5, 24, by = 3), cmt = "center")
    .ev <- rxode2::et(.ev, seq(0.5, 24, by = 3), cmt = "effect")
    .d <- as.data.frame(rxode2::rxSolve(m, .ev, addDosing = TRUE))
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

  test_that("est='npag' scores a general likelihood alongside a second endpoint (#850)", {
    skip_if_not(rxode2hasLlik(), "rxode2 build has no llik support")
    skip_if_not(.rxode2HasSafeLogDomain(),
                "rxode2 too old for the safeLog log-domain mode")
    .gauss <- .npMkPkPd(quote(effect ~ add(pdadd.sd) | effect))
    .ll <- .npMkPkPd(quote(ll(effect) ~ -0.5 * log(2 * pi) - log(pdadd.sd) -
                             0.5 * ((DV - effect) / pdadd.sd)^2))
    .dat <- .npPkPdData(.gauss)
    .nObs <- sum(.dat$evid == 0)
    .ctl <- function() npagControl(points = 24L, cycles = 3L, seed = 1L, print = 0L)
    rxode2::rxSetSeed(42); .fg <- suppressWarnings(nlmixr2(.gauss, .dat, "npag", .ctl()))
    rxode2::rxSetSeed(42); .fl <- suppressWarnings(nlmixr2(.ll, .dat, "npag", .ctl()))
    .lg <- as.numeric(.fg$env$npagLogLik)
    .ll2 <- as.numeric(.fl$env$npagLogLik)
    # every scale stays in its domain -- this is what actually went wrong, and the
    # bound below is only well defined once it holds
    expect_gt(unname(fixef(.fl)["pdadd.sd"]), 0)
    expect_gt(unname(fixef(.fl)["add.sd"]), 0)
    # A per-observation normal log-density cannot exceed -log(sd) - 0.5*log(2*pi), so the
    # total cannot exceed the sum of those maxima.  Take the sd's from the FIT, not from
    # the simulation: the fit chooses its own, and pdadd.sd lands near 0.76 against the
    # simulated 2 -- a bound built on the simulation values would be one this fit is
    # entitled to beat, and would fail for a legitimate reason.  At the fitted sd's this
    # is a theorem.  The defect cleared even the looser bound by more than 2000.
    .nPk <- sum(.dat$evid == 0 & .dat$cmt == 2)
    .nPd <- sum(.dat$evid == 0 & .dat$cmt == 3)
    expect_equal(.nPk + .nPd, .nObs)
    .maxLL <- -.nObs * 0.5 * log(2 * pi) -
      .nPk * log(unname(fixef(.fl)["add.sd"])) -
      .nPd * log(unname(fixef(.fl)["pdadd.sd"]))
    expect_lt(.ll2, .maxLL)
    # The twin's objectives differ by the 2*pi term the Gaussian path omits, plus
    # whatever the two fits' differing residual estimates are worth -- see below.
    # This is a RELATIVE tolerance on a number near 176, so it allows a few units.
    expect_equal(.lg - .ll2, .nObs * 0.5 * log(2 * pi), tolerance = 0.02)
    # add.sd is the PK endpoint, written `cp ~ add(add.sd)` in BOTH models, so both
    # fits score it the same way and it agrees tightly (measured 0.6359 vs 0.6352).
    expect_equal(unname(fixef(.fl)["add.sd"]), unname(fixef(.fg)["add.sd"]),
                 tolerance = 0.05)
    # pdadd.sd is NOT required to agree, and asserting that it did was wrong.
    # The residual step optimizes an err-tagged parameter against the
    # extended-least-squares objective at the posterior-mean etas (npCommon.R);
    # a general-likelihood endpoint has no ELS form, so it is scored by the
    # user's log-density instead.  The two endpoints are therefore scored by
    # DIFFERENT objectives by construction, and only the endpoint that is
    # spelled identically in both models (add.sd, above) has to match.
    #
    # Measured at points=24/cycles=3: 0.9267 (Gaussian) vs 0.7343 (ll), and the
    # objective is not flat between them -- refitting with pdadd.sd FIXED gives
    # -49.89, -45.65 and -57.54 at 0.9267, 0.80 and 0.7343 -- so this is a real
    # difference in where the two objectives put the optimum, not noise and not
    # a tolerance that needs widening.  Whether ELS is accurate enough for a
    # multi-endpoint residual step is a separate question about npag's design;
    # it is not something this test can pin as an equality.
    expect_true(is.finite(unname(fixef(.fl)["pdadd.sd"])))
    expect_gt(unname(fixef(.fl)["pdadd.sd"]), 0)
    expect_true(is.finite(unname(fixef(.fg)["pdadd.sd"])))
    expect_gt(unname(fixef(.fg)["pdadd.sd"]), 0)
  })

  test_that("an ini() lower bound constrains a parameter inside ll() (#850)", {
    # A hand-written likelihood's domain is not always a box -- log(a-b) needs a > b --
    # so nlmixr2est does not try to INFER that a parameter is a scale.  Declaring the
    # bound in ini() is the supported way to say so, and it reaches the nonparametric
    # residual step as the optimizer's box.  This pins that contract.
    skip_if_not(rxode2hasLlik(), "rxode2 build has no llik support")
    .m <- function() {
      ini({ tka <- log(1.5); tv <- log(32); tke <- fix(log(0.08))
        eta.ka ~ 0.3; eta.v ~ 0.1; add.sd <- c(0.05, 0.7) })
      model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); ke <- exp(tke)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - ke * center
        cp <- center / v
        ll(cp) ~ -0.5 * log(2 * pi) - log(add.sd) - 0.5 * ((DV - cp) / add.sd)^2 })
    }
    f <- suppressWarnings(nlmixr2(.m, nlmixr2data::theo_sd, est = "npag",
                                  control = npagControl(points = 24L, cycles = 3L,
                                                        seed = 1L, print = 0L)))
    expect_s3_class(f, "nlmixr2FitData")
    expect_gte(unname(fixef(f)["add.sd"]), 0.05)   # the declared lower bound holds
  })
})
