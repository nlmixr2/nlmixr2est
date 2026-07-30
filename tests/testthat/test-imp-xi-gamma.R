# Tests for the NONMEM-style per-subject xi statistic and the individual-gamma
# controller built on it.
#
# xi (NM7 Technical Guide eq. 1.73/1.75) is the mean importance weight normalized
# so the proposal kernel is 1 at its center:
#
#   xi_i = (1/r) sum_k [pi(eta_k) / pi(center)] / e(eta_k)
#
# It equals 1 when the proposal matches the posterior shape, is < 1 when the
# proposal is over-dispersed, and is > 1 when the posterior has heavier tails
# than the Gaussian proposal.  It is a DIFFERENT statistic from the Kish
# effective sample size (impNeffFrac) and the two are not comparable.
nmTest({

  .xiModel <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.cl ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("xi: per-subject xi and its trace are exposed and well-formed", {
    .d <- nlmixr2data::theo_sd
    .f <- suppressWarnings(nlmixr2(.xiModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 5L,
                                                 isample = 300L, covMethod = "")))
    .E <- .f$env
    # one xi per subject, all finite and strictly positive
    expect_false(is.null(.E$impXi))
    expect_length(.E$impXi, length(unique(.d$ID)))
    expect_true(all(is.finite(.E$impXi)))
    expect_true(all(.E$impXi > 0))
    # the trace has one entry per EM iteration actually run, alongside the
    # Kish-ESS trace it must not be confused with
    expect_length(.E$impXiTrace, .E$impIter)
    expect_length(.E$impNeffFrac, .E$impIter)
    expect_true(all(is.finite(.E$impXiTrace)))
    # xi and the Kish ESS fraction are different statistics: on a real fit they
    # must not coincide, which is exactly why $runInfo has to name which is which
    expect_false(isTRUE(all.equal(.E$impXiTrace, .E$impNeffFrac)))
  })

  test_that("xi decreases monotonically as the proposal is widened (gamma)", {
    # This is the property that makes xi targetable by the individual-gamma
    # controller: widening the proposal must lower xi, monotonically.  Compared
    # at iteration 1 across fits so every run sees identical starting parameters.
    .d <- nlmixr2data::theo_sd
    .xi1 <- vapply(c(1.0, 2.0, 4.0), function(g) {
      .f <- suppressWarnings(nlmixr2(.xiModel, .d, "impmap",
                                     impmapControl(print = 0L, nIter = 1L,
                                                   isample = 500L, gamma = g,
                                                   covMethod = "")))
      .f$env$impXiTrace[1]
    }, numeric(1))
    expect_true(all(is.finite(.xi1)))
    expect_true(all(diff(.xi1) < 0))
  })

  test_that("gammaMethod: control accepts both modes and defaults to global", {
    expect_equal(impmapControl()$gammaMethod, "global")
    expect_equal(impmapControl(gammaMethod = "individual")$gammaMethod, "individual")
    expect_error(impmapControl(gammaMethod = "nonsense"))
    # round-trips through the control constructor (do.call re-entry)
    .c <- impmapControl(gammaMethod = "individual")
    expect_equal(do.call(impmapControl, .c)$gammaMethod, "individual")
    # stripped when down-converting to a plain foceiControl
    expect_true("gammaMethod" %in% .impmapIsControlNames)
  })

  test_that("gammaMethod='global' keeps one shared scale", {
    .d <- nlmixr2data::theo_sd
    .f <- suppressWarnings(nlmixr2(.xiModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 10L,
                                                 covMethod = "",
                                                 gammaMethod = "global")))
    .E <- .f$env
    expect_equal(.E$impGammaMethod, "global")
    # every subject on the same scale, equal to the reported scalar
    expect_equal(length(unique(round(.E$impGammaInd, 12))), 1L)
    expect_equal(unname(.E$impGammaInd[1]), unname(.E$impGammaUsed), tolerance = 1e-10)
  })

  test_that("gammaMethod='individual' gives per-subject scales that target iaccept", {
    .d <- nlmixr2data::theo_sd
    .f <- suppressWarnings(nlmixr2(.xiModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 30L,
                                                 covMethod = "",
                                                 gammaMethod = "individual")))
    .E <- .f$env
    expect_equal(.E$impGammaMethod, "individual")
    # subjects genuinely diverge -- this is the whole point of the mode
    expect_gt(length(unique(round(.E$impGammaInd, 8))), 1L)
    # every scale respects the ISCALE bounds
    expect_true(all(.E$impGammaInd >= 0.1 - 1e-8))
    expect_true(all(.E$impGammaInd <= 10 + 1e-8))
    expect_true(all(is.finite(.E$impGammaInd)))
    # the controller drives mean xi onto the iaccept target (default 0.4).  The
    # global rule leaves it near 1 on this near-Gaussian model, so this is a
    # genuine discriminator between the two laws, not a tautology.
    expect_equal(tail(.E$impXiTrace, 1), 0.4, tolerance = 0.1)
  })

  test_that("gammaMethod does not move the estimates, only the variance", {
    # The importance weights correct for gamma, so both laws must agree on the
    # parameters even though they sample at very different proposal widths.
    .d <- nlmixr2data::theo_sd
    .g <- suppressWarnings(nlmixr2(.xiModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 25L, covMethod = "",
                                                 gammaMethod = "global")))
    .i <- suppressWarnings(nlmixr2(.xiModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 25L, covMethod = "",
                                                 gammaMethod = "individual")))
    # the two ran at materially different scales ...
    expect_gt(.i$env$impGammaUsed, 1.2 * .g$env$impGammaUsed)
    # ... but agree on the fixed effects and the objective
    expect_equal(unname(.i$theta), unname(.g$theta), tolerance = 0.02)
    expect_equal(.i$env$impObj, .g$env$impObj, tolerance = 0.5)
  })

  test_that("individual gamma respects a tightened iscaleMax clamp", {
    .d <- nlmixr2data::theo_sd
    .f <- suppressWarnings(nlmixr2(.xiModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 20L, covMethod = "",
                                                 gammaMethod = "individual",
                                                 iscaleMax = 1.2)))
    # the unclamped controller wants gamma ~ 1.8 here, so 1.2 must bind
    expect_true(all(.f$env$impGammaInd <= 1.2 + 1e-8))
    expect_gt(max(.f$env$impGammaInd), 1.0)
  })

  test_that("covMethod='imp' uses the converged scales under individual gamma", {
    # impComputeCov used to read impGammaProp() -- the control's INITIAL gamma --
    # so with gammaMethod="individual" it evaluated the covariance with a single
    # scale of 1.0 while the fit had converged near 1.8 per subject.  The
    # importance weights correct for gamma, so that was consistent rather than
    # biased, but it sampled from a proposal the fit had already rejected as
    # poorly matched, which inflates Monte-Carlo noise in the FD Hessian.  This
    # pins that the covariance is built and stays well-formed on that path.
    .d <- nlmixr2data::theo_sd
    .f <- suppressWarnings(nlmixr2(.xiModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 15L,
                                                 covMethod = "imp",
                                                 gammaMethod = "individual")))
    .E <- .f$env
    # the scales really did move away from the control's initial gamma = 1
    expect_gt(min(.E$impGammaInd), 1.05)
    .cv <- as.matrix(.E$cov)
    expect_true(all(is.finite(.cv)))
    expect_equal(.cv, t(.cv), tolerance = 1e-8)          # symmetric
    expect_true(all(diag(.cv) > 0))                       # positive diagonal
    expect_true(all(is.finite(.E$impSeTheta)))
    expect_true(all(.E$impSeTheta > 0))
  })

  test_that("individual gamma is stable at high eta dimension (no limit cycle)", {
    # The damping exponent has to scale with neta.  For a Gaussian posterior
    # xi = gamma^(-neta/2), so gamma' = gamma*(xi/iaccept)^p linearizes to an
    # error multiplier lambda = 1 - p*neta/2.  A FIXED p = 1/2 gives
    # lambda = -1 at neta = 8, i.e. a period-2 limit cycle that never settles
    # (observed: gamma 1.32/1.24/1.30/1.23..., xi 0.36/0.45/0.37/0.46...).
    # p = 2/neta gives lambda = 0 for every dimension.  This test fails on the
    # fixed-exponent controller and passes on the dimension-aware one.
    skip_on_cran()
    .testSeed(7); rxode2::rxSetSeed(7)
    .mk <- function() {
      tt <- c(0.25, 1, 2, 4, 8, 16, 24)
      do.call(rbind, lapply(1:12, function(id) {
        ka <- exp(0.4 + stats::rnorm(1, 0, .3)); cl <- exp(1 + stats::rnorm(1, 0, .3))
        v <- exp(3.4 + stats::rnorm(1, 0, .3))
        cp <- 100 * ka / (v * (ka - cl / v)) * (exp(-cl / v * tt) - exp(-ka * tt))
        cp <- pmax(cp, 1e-3) * exp(stats::rnorm(length(tt), 0, .1))
        rbind(data.frame(id = id, time = 0, dv = NA_real_, amt = 100, evid = 1, cmt = "depot"),
              data.frame(id = id, time = tt, dv = cp, amt = 0, evid = 0, cmt = "cen"))
      }))
    }
    .d <- .mk(); .d <- .d[order(.d$id, .d$time, -.d$evid), ]
    m8 <- function() {
      ini({
        tka <- 0.4; tcl <- 1; tv <- 3.4; tq <- 0.8; tv2 <- 3; tf <- 0; tlag <- -2; tke0 <- -1
        eta.ka ~ .3; eta.cl ~ .3; eta.v ~ .3; eta.q ~ .3
        eta.v2 ~ .3; eta.f ~ .2; eta.lag ~ .2; eta.ke0 ~ .2
        add.sd <- 0.3
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        q <- exp(tq + eta.q); v2 <- exp(tv2 + eta.v2); f <- exp(tf + eta.f)
        lag <- exp(tlag + eta.lag); ke0 <- exp(tke0 + eta.ke0)
        d/dt(depot) <- -ka * depot
        d/dt(cen) <- ka * depot - cl / v * cen - q / v * cen + q / v2 * peri + 0 * lag + 0 * ke0
        d/dt(peri) <- q / v * cen - q / v2 * peri
        cp <- f * cen / v
        cp ~ add(add.sd)
      })
    }
    .f <- suppressWarnings(nlmixr2(m8, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 18L, isample = 200L,
                                                 covMethod = "", gammaMethod = "individual")))
    expect_equal(nrow(.f$omega), 8L)
    .tail <- tail(.f$env$impXiTrace, 10)
    # a period-2 limit cycle has lag-1 autocorrelation ~ -1; a settled trace
    # with Monte-Carlo noise sits well above that
    expect_gt(stats::cor(head(.tail, -1), .tail[-1]), -0.8)
    # and the swing must be small relative to the target it is holding
    expect_lt(diff(range(.tail)), 0.15)
  })

  test_that("xi is near 1 for a well-matched proposal on a near-Gaussian model", {
    # gamma = 1 makes the proposal the Laplace approximation itself, so on a
    # model whose individual posterior is close to Gaussian xi should sit near 1.
    # Deliberately a loose band -- this pins the scale/orientation of the
    # statistic (a sign error or a missing normalizer would blow it out by
    # orders of magnitude), not its precise value.
    .d <- nlmixr2data::theo_sd
    .f <- suppressWarnings(nlmixr2(.xiModel, .d, "impmap",
                                   impmapControl(print = 0L, nIter = 1L,
                                                 isample = 500L, gamma = 1.0,
                                                 covMethod = "")))
    expect_true(.f$env$impXiTrace[1] > 0.5)
    expect_true(.f$env$impXiTrace[1] < 2.0)
  })

})
