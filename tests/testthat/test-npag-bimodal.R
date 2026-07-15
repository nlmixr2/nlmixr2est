## Bimodal recovery -- the defining nonparametric capability: NPAG recovers a
## two-subpopulation (fast/slow absorption) parameter distribution that a
## parametric normal random-effect model would collapse into a single mode.
## Simulated data with a known bimodal Ka; check the support points cluster into
## two groups near the true values.  Real fit -> weekly slow batch.

nmTest({
  test_that("est='npag' recovers a bimodal Ka distribution", {
    skip_if_not_installed("rxode2")
    set.seed(42)
    nEach <- 15L
    kaSlow <- 0.7; kaFast <- 3.5           # two subpopulations
    vTrue <- 30; keTrue <- 0.1
    kaTrue <- c(rep(kaSlow, nEach), rep(kaFast, nEach))
    N <- length(kaTrue)

    simMod <- rxode2::rxode2({
      d/dt(depot) <- -KA * depot
      d/dt(center) <- KA * depot - KE * center
      cp <- center / V
    })
    ev <- rxode2::et(amt = 100, cmt = "depot")
    for (.t in c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 24)) ev <- rxode2::et(ev, .t)
    pars <- data.frame(KA = kaTrue, V = vTrue, KE = keTrue)
    s <- rxode2::rxSolve(simMod, pars, ev, returnType = "data.frame")
    # the subject column may be named "id" or "sim.id" depending on rxode2 version
    .idCol <- names(s)[grepl("id$", names(s), ignore.case = TRUE)][1]
    s$.id <- s[[.idCol]]
    # assemble a NONMEM-style dataset with additive residual error
    obs <- s[s$time > 0, ]
    obs$DV <- obs$cp + rnorm(nrow(obs), 0, 0.3)
    dat <- data.frame(ID = obs$.id, TIME = obs$time, DV = obs$DV,
                      AMT = 0, EVID = 0, CMT = 2)
    dose <- data.frame(ID = unique(obs$.id), TIME = 0, DV = 0,
                       AMT = 100, EVID = 1, CMT = 1)
    dat <- rbind(dose, dat)
    dat <- dat[order(dat$ID, dat$TIME, -dat$EVID), ]

    npMod <- function() {
      ini({ tka <- log(1.5); tv <- log(30); tke <- log(0.1)
        eta.ka ~ 0.5; eta.v ~ 0.05; eta.ke ~ 0.05; add.sd <- 0.3 })
      model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); ke <- exp(tke + eta.ke)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - ke * center
        cp <- center / v
        cp ~ add(add.sd) })
    }
    f <- nlmixr2(npMod, dat, est = "npag",
                 control = npagControl(points = 512L, cycles = 40L, gammaOptimize = FALSE))
    expect_s3_class(f, "nlmixr2FitData")

    # support-point Ka = exp(tka + eta.ka); the mixture should have mass at BOTH
    # subpopulations, not a single averaged mode.
    tka <- as.numeric(f$theta[["tka"]])
    kaSupport <- exp(tka + f$env$npagSupport[, 1])
    wt <- f$env$npagWeights
    mid <- sqrt(kaSlow * kaFast)                      # geometric-mean split point
    wSlow <- sum(wt[kaSupport < mid])
    wFast <- sum(wt[kaSupport >= mid])
    # both modes carry substantial weight (a unimodal fit would starve one side)
    expect_gt(wSlow, 0.2)
    expect_gt(wFast, 0.2)
    # the weighted cluster means land near the true Ka values
    mSlow <- sum(wt[kaSupport < mid] * kaSupport[kaSupport < mid]) / wSlow
    mFast <- sum(wt[kaSupport >= mid] * kaSupport[kaSupport >= mid]) / wFast
    expect_lt(abs(mSlow - kaSlow), 0.5)
    expect_lt(abs(mFast - kaFast), 1.2)
    # the two recovered modes are clearly separated
    expect_gt(mFast / mSlow, 2)
  })
})
