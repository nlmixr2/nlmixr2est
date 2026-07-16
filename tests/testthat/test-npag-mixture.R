## Mixture (sub-population) support for NPAG -- a mix() model splits each subject
## into per-component pseudo-subjects.  NPAG marginalizes the conditional
## likelihood over the components using the mixture proportions
## (p(y_i | phi) = sum_m mixProb_m * p(y_i | phi, component m)); when a proportion
## is estimated it is updated by an in-cycle EM step (support points + weights
## held fixed), and when fixed it is held at its ini value.  Real fit -> weekly
## slow batch.

nmTest({
  test_that("est='npag' holds a FIXED mix() proportion and marginalizes over it", {
    # p1 fixed => no EM update; the two fixed values give different marginal
    # likelihoods (a proportion-independent objf would mean the components were
    # never combined -- the pre-marginalization bug).
    .mix5 <- function() {
      ini({ tka <- log(1.5); tv <- log(32); tcl1 <- log(2); tcl2 <- log(0.5)
        eta.ka ~ 0.3; eta.v ~ 0.1; p1 <- fix(0.5); add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v)
        cl <- mix(exp(tcl1), p1, exp(tcl2)); ke <- cl / v
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - ke * center
        cp <- center / v; cp ~ add(add.sd) })
    }
    .mix9 <- function() {
      ini({ tka <- log(1.5); tv <- log(32); tcl1 <- log(2); tcl2 <- log(0.5)
        eta.ka ~ 0.3; eta.v ~ 0.1; p1 <- fix(0.9); add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v)
        cl <- mix(exp(tcl1), p1, exp(tcl2)); ke <- cl / v
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - ke * center
        cp <- center / v; cp ~ add(add.sd) })
    }
    .ctl <- npagControl(points = 32L, cycles = 3L, seed = 1L, gammaOptimize = FALSE, muExpand = FALSE)
    f5 <- nlmixr2(.mix5, nlmixr2data::theo_sd, est = "npag", control = .ctl)
    f9 <- nlmixr2(.mix9, nlmixr2data::theo_sd, est = "npag", control = .ctl)
    expect_s3_class(f5, "nlmixr2FitData")
    expect_equal(as.numeric(f5$theta[["p1"]]), 0.5)   # fixed -> held
    expect_equal(as.numeric(f9$theta[["p1"]]), 0.9)
    expect_false(isTRUE(all.equal(as.numeric(f5$objf), as.numeric(f9$objf))))
  })

  test_that("est='npb' fits a mix() model and samples the mixture proportion", {
    .mixMod <- function() {
      ini({ tka <- log(1.5); tv <- log(32); tcl1 <- log(2); tcl2 <- log(0.5)
        eta.v ~ 0.1; p1 <- 0.5; add.sd <- 0.7 })
      model({ ka <- exp(tka); v <- exp(tv + eta.v)
        cl <- mix(exp(tcl1), p1, exp(tcl2)); ke <- cl / v
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - ke * center
        cp <- center / v; cp ~ add(add.sd) })
    }
    f <- nlmixr2(.mixMod, nlmixr2data::theo_sd, est = "npb",
                 control = npbControl(points = 20L, burnin = 20L, nsamp = 20L, seed = 1L))
    expect_s3_class(f, "nlmixr2FitData")
    expect_true(is.finite(as.numeric(f$objf)))
    # the mixture proportions are sampled (Dirichlet Gibbs) -> a valid probability
    # vector, reported in $env$npbMixProb, that has moved off the 0.5 start.
    .p <- f$env$npbMixProb
    expect_equal(length(.p), 2L)
    expect_equal(sum(.p), 1, tolerance = 1e-6)
    expect_true(all(.p > 0 & .p < 1))
    expect_true(as.numeric(f$theta[["p1"]]) > 0 && as.numeric(f$theta[["p1"]]) < 1)
  })

  test_that("est='npag' estimates the mix() proportion via the in-cycle EM update", {
    skip_if_not_installed("rxode2")
    set.seed(7)
    N <- 40L
    comp <- rbinom(N, 1, 0.3)                  # 1 = fast-clearance subpopulation
    clTrue <- ifelse(comp == 1, 2, 0.5)
    sim <- rxode2::rxode2({
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - (cl / v) * center
      cp <- center / v
    })
    ev <- rxode2::et(amt = 100, cmt = "depot")
    for (.t in c(0.5, 1, 2, 4, 6, 8, 12, 24)) ev <- rxode2::et(ev, .t)
    s <- rxode2::rxSolve(sim, data.frame(ka = 1.2, v = 32, cl = clTrue), ev,
                         returnType = "data.frame")
    .idc <- names(s)[grepl("id$", names(s), ignore.case = TRUE)][1]
    s$.id <- s[[.idc]]; obs <- s[s$time > 0, ]
    obs$DV <- obs$cp + rnorm(nrow(obs), 0, 0.3)
    dat <- rbind(
      data.frame(ID = unique(obs$.id), TIME = 0, DV = 0, AMT = 100, EVID = 1, CMT = 1),
      data.frame(ID = obs$.id, TIME = obs$time, DV = obs$DV, AMT = 0, EVID = 0, CMT = 2))
    dat <- dat[order(dat$ID, dat$TIME, -dat$EVID), ]

    mixMod <- function() {
      ini({ tka <- log(1.2); tv <- log(32); tcl1 <- log(0.5); tcl2 <- log(2)
        eta.v ~ 0.05; p1 <- 0.5; add.sd <- 0.3 })
      model({ ka <- exp(tka); v <- exp(tv + eta.v)
        cl <- mix(exp(tcl1), p1, exp(tcl2)); ke <- cl / v
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - ke * center
        cp <- center / v; cp ~ add(add.sd) })
    }
    # muExpand=FALSE: this test validates the proportion EM with a fixed structural
    # model (tka/tcl held); mu-expanding them adds free parameters that weaken the
    # proportion identification.  Mixture + mu-expansion is covered separately.
    f <- nlmixr2(mixMod, dat, est = "npag",
                 control = npagControl(points = 64L, cycles = 20L, seed = 1L,
                                       gammaOptimize = FALSE, muExpand = FALSE))
    expect_s3_class(f, "nlmixr2FitData")
    # the slow-clearance proportion (p1) recovers the simulated fraction (0.70)
    # from a deliberately-wrong 0.50 start -- the EM update moved it.
    expect_equal(as.numeric(f$theta[["p1"]]), mean(comp == 0), tolerance = 0.1)
  })

  test_that("est='npag' mu-expands a mixture model's structural thetas (multi-eta)", {
    # a mixture model with non-mu structural thetas (tka, and the component values
    # tcl1/tcl2) -- with muExpand each gets a fixed-omega pseudo-eta.  This exercises
    # >1 eta in a mixture, which previously errored "initial 'omega' matrix inverse
    # is non-positive definite"; the impSetOmega PD guard fixes it.  The mixture
    # proportion (p1) is NOT injected -- it is estimated by the in-cycle EM.
    .m <- function() {
      ini({ tka <- log(1.5); tv <- log(32); tcl1 <- log(2); tcl2 <- log(0.5)
        eta.v ~ 0.1; p1 <- 0.5; add.sd <- 0.7 })
      model({ ka <- exp(tka); v <- exp(tv + eta.v)
        cl <- mix(exp(tcl1), p1, exp(tcl2)); ke <- cl / v
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - ke * center
        cp <- center / v; cp ~ add(add.sd) })
    }
    f <- nlmixr2(.m, nlmixr2data::theo_sd, est = "npag",
                 control = npagControl(points = 32L, cycles = 3L, seed = 1L,
                                       gammaOptimize = FALSE, calcTables = FALSE,
                                       muExpand = TRUE))
    expect_true(is.finite(as.numeric(f$objf)))
    # tka and the component clearances gained injected (collapsed) etas ...
    expect_true(all(c("eta.tka", "eta.tcl1", "eta.tcl2") %in% rownames(f$omega)))
    expect_lt(f$omega["eta.tka", "eta.tka"], 1e-4)         # collapsed to a fixed effect
    # ... while the proportion is not injected and is a valid probability
    expect_false("eta.p1" %in% rownames(f$omega))
    expect_true(as.numeric(f$theta[["p1"]]) > 0 && as.numeric(f$theta[["p1"]]) < 1)
  })
})
