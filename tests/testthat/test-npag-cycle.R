## NPAG adaptive-grid cycle (npagCycle_) end-to-end on Theophylline: Sobol grid
## -> Psi -> Burke IPM -> condensation -> expansion -> convergence.  Runs a real
## (if modest) fit, so it lives in a weekly slow batch.

nmTest({
  test_that("npagCycle_ converges to a sane discrete distribution on theo_sd", {
    npMod <- function() {
      ini({ tka <- log(1.5); tv <- log(32); tke <- log(0.08)
        eta.ka ~ 0.3; eta.v ~ 0.1; eta.ke ~ 0.1
        add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); ke <- exp(tke + eta.ke)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - ke * center
        cp <- center / v
        cp ~ add(add.sd) })
    }
    ui <- rxode2::assertRxUi(npMod)
    ctl <- npagControl()
    dat <- nlmixr2data::theo_sd
    N <- length(unique(dat$ID))

    .npInnerSetup(ui, dat, matrix(0, N, 3L), ctl)
    on.exit(.npInnerFree(), add = TRUE)
    box <- .npEtaBox(ui, ctl)
    cores <- as.integer(rxode2::getRxThreads())

    r <- npagCycle_(box$lower, box$upper, points = 256L, cycles = 20L, cores = cores)

    # weights form a probability vector
    expect_equal(sum(r$weights), 1, tolerance = 1e-6)
    expect_true(all(r$weights >= -1e-10))
    # Caratheodory: the NPML has at most N support points
    nspp <- nrow(r$support)
    expect_true(nspp >= 1L && nspp <= N)
    expect_equal(ncol(r$support), 3L)
    # support points stay within the box
    for (j in 1:3) {
      expect_true(all(r$support[, j] >= box$lower[j] - 1e-8 &
                        r$support[, j] <= box$upper[j] + 1e-8))
    }
    # objective / -2LL are finite and in a sane range for this 12-subject PK model
    expect_true(is.finite(r$objf))
    m2ll <- -2 * r$objf
    expect_true(m2ll > 50 && m2ll < 5000)
    expect_true(r$converged || r$cycles == 20L)

    # the weighted-mean eta is a small deviation (etas are individual offsets)
    wmean <- colSums(r$support * r$weights)
    expect_true(all(abs(wmean) < 2))
  })

  test_that("gamma optimization improves the objective and returns a finite gamma", {
    npMod <- function() {
      ini({ tka <- log(1.5); tv <- log(32); tke <- log(0.08)
        eta.ka ~ 0.3; eta.v ~ 0.1; eta.ke ~ 0.1
        add.sd <- 1.4 })   # deliberately too large -> gamma should shrink it (< 1)
      model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); ke <- exp(tke + eta.ke)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - ke * center
        cp <- center / v
        cp ~ add(add.sd) })
    }
    ui <- rxode2::assertRxUi(npMod)
    ctl <- npagControl()
    dat <- nlmixr2data::theo_sd
    N <- length(unique(dat$ID))
    .npInnerSetup(ui, dat, matrix(0, N, 3L), ctl)
    on.exit(.npInnerFree(), add = TRUE)
    box <- .npEtaBox(ui, ctl)
    cores <- as.integer(rxode2::getRxThreads())

    r0 <- npagCycle_(box$lower, box$upper, points = 256L, cycles = 20L,
                     cores = cores, gammaOptimize = FALSE)
    rg <- npagCycle_(box$lower, box$upper, points = 256L, cycles = 20L,
                     cores = cores, gammaOptimize = TRUE)

    expect_true(is.finite(rg$gamma) && rg$gamma > 0)
    expect_equal(sum(rg$weights), 1, tolerance = 1e-6)
    # optimizing gamma cannot do worse than holding it at 1
    expect_gte(rg$objf, r0$objf - 1e-6)
    # add.sd = 1.4 is inflated for theo, so the optimum gamma multiplier < 1
    expect_lt(rg$gamma, 1)
    # the objf(gamma) landscape has an interior peak at the optimized gamma, i.e.
    # the log-variance penalty scales with gamma -- guards the "larger gamma always
    # wins" regression.  Checked in the found gamma's own neighborhood since the
    # peak location depends on the (inflated) starting add.sd.
    g <- rg$support; gOpt <- rg$gamma
    fOpt <- npObjAtGamma_(g, cores, gOpt)
    expect_gt(fOpt, npObjAtGamma_(g, cores, gOpt * 0.2))
    expect_gt(fOpt, npObjAtGamma_(g, cores, gOpt * 5.0))
  })

  test_that("censored (BLQ) observations flow through the gamma path", {
    npMod <- function() {
      ini({ tka <- log(1.5); tv <- log(32); tke <- log(0.08)
        eta.ka ~ 0.3; eta.v ~ 0.1; eta.ke ~ 0.1; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); ke <- exp(tke + eta.ke)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - ke * center
        cp <- center / v
        cp ~ add(add.sd) })
    }
    ui <- rxode2::assertRxUi(npMod)
    ctl <- npagControl()
    # censor low observations as BLQ (CENS = 1, DV = LOQ) -- M3 likelihood
    dat <- nlmixr2data::theo_sd
    LOQ <- 2.5
    dat$CENS <- 0L
    .obs <- dat$EVID == 0 & dat$DV < LOQ
    dat$CENS[.obs] <- 1L; dat$DV[.obs] <- LOQ
    expect_gt(sum(.obs), 10L)     # a meaningful number of BLQ points
    N <- length(unique(dat$ID))
    cores <- as.integer(rxode2::getRxThreads())
    .npInnerSetup(ui, dat, matrix(0, N, 3L), ctl)
    on.exit(.npInnerFree(), add = TRUE)
    box <- .npEtaBox(ui, ctl)

    rg <- npagCycle_(box$lower, box$upper, points = 256L, cycles = 20L,
                     cores = cores, gammaOptimize = TRUE)
    expect_equal(sum(rg$weights), 1, tolerance = 1e-6)
    expect_true(is.finite(rg$objf) && rg$gamma > 0 && rg$gamma < 5)
    expect_true(nrow(rg$support) >= 1L && nrow(rg$support) <= N)
    # the censored likelihood scales with gamma (npResidScale reaches
    # doCensNormal1): objf(gamma) peaks in the interior, not at the boundary.
    g <- rg$support
    fPk <- npObjAtGamma_(g, cores, 1.0)
    expect_gt(fPk, npObjAtGamma_(g, cores, 0.5))
    expect_gt(fPk, npObjAtGamma_(g, cores, 2.5))
  })
})
