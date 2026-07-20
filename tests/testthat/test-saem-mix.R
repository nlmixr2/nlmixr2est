nmTest({
  test_that(".configsaem only floors the i0 (no-eta) phi variance at 1.0 for mixture fits (nMix>1)", {
    # i0 indexes phi positions with NO associated eta at all (ordinary
    # fixed-effect/population-only thetas, e.g. tka below) -- not
    # "non-mu-referenced etas". Capture .configsaem()'s returned cfg (which
    # includes minv and i0 directly) via trace()'s exit hook, aborting right
    # after .configsaem returns so this test doesn't pay for a full
    # stochastic SAEM fit. The tracer runs inside .configsaem's own call
    # frame, so it writes to .GlobalEnv (always lexically reachable)
    # rather than a test-local environment.
    withr::defer(if (exists(".sddCfgCapture", envir = .GlobalEnv)) {
      rm(".sddCfgCapture", envir = .GlobalEnv)
    })
    .captureCfg <- function() {
      trace(".configsaem",
            exit = quote({
              assign(".sddCfgCapture", returnValue(), envir = .GlobalEnv)
              stop("test-capture-exit")
            }),
            print = FALSE, where = asNamespace("nlmixr2est"))
      withr::defer(suppressMessages(untrace(".configsaem", where = asNamespace("nlmixr2est"))),
                   envir = parent.frame())
    }

    d <- do.call(rbind, lapply(1:8, function(i) {
      times <- c(0.5, 1, 2, 4, 8)
      data.frame(ID = i, TIME = c(0, times), AMT = c(100, rep(0, length(times))),
                 EVID = c(1, rep(0, length(times))), DV = c(0, rep(1, length(times))),
                 CMT = c(1, rep(2, length(times))))
    }))

    # Non-mixture model: tka has no eta at all (i0 element).
    one.compartment.i0 <- function() {
      ini({
        tka <- log(1.5)
        tcl <- log(2.0)
        tv <- log(20)
        eta.cl ~ 0.01
        eta.ka2 ~ 0.09
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka) + eta.ka2
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .captureCfg()
    suppressWarnings(suppressMessages(try(
      .nlmixr(one.compartment.i0, d, est = "saem",
              saemControl(print = 0, nBurn = 1, nEm = 1, calcTables = FALSE)),
      silent = TRUE
    )))
    expect_true(length(.sddCfgCapture$i0) >= 1)
    expect_equal(unname(.sddCfgCapture$minv[.sddCfgCapture$i0 + 1L]),
                 rep(1e-20, length(.sddCfgCapture$i0)))
    rm(".sddCfgCapture", envir = .GlobalEnv)

    # Mixture model (nMix=2): same no-eta tka, but mixProb has length > 1.
    one.compartment.mix.i0 <- function() {
      ini({
        tka <- log(1.5)
        tcl1 <- log(1.0)
        tcl2 <- log(5.0)
        tv <- log(20)
        p1 <- 0.5
        eta.cl ~ 0.01
        eta.v ~ 0.01
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .captureCfg()
    suppressWarnings(suppressMessages(try(
      .nlmixr(one.compartment.mix.i0, d, est = "saem",
              saemControl(print = 0, nBurn = 1, nEm = 1, calcTables = FALSE)),
      silent = TRUE
    )))
    expect_true(length(.sddCfgCapture$i0) >= 1)
    expect_equal(unname(.sddCfgCapture$minv[.sddCfgCapture$i0 + 1L]),
                 rep(1.0, length(.sddCfgCapture$i0)))
  })

  test_that("SAEM warns when a mixture probability estimate collapses/needs rescaling", {
    withr::local_options(list(warn = 0))
    # Directly exercise the clamp/warn logic via the internal helper that
    # runs at the end of every SAEM mixture fit (.getSaemTheta), rather than
    # forcing a real SAEM fit to a collapsed solution (slow/flaky). The `ui`
    # objects below are real rxode2 UI objects (not hand-rolled stand-ins),
    # so this exercises the actual production code path; only the raw SAEM
    # optimizer output (`env$saem`) is mocked.
    #
    # Internally, mixProb components each individually stay in [0, 1] (they
    # are updated by a convex-combination stochastic-approximation step), so
    # a single component near 0/1 (e.g. 1e-8) is only ~1e-6 away from the
    # clamp boundary -- too small to cross the "real change" threshold below.
    # The scenario that plausibly produces a large, warning-worthy change is
    # multiple *independently* estimated probabilities whose sum drifts to
    # >= 1 (non-identifiability across >2 components), which is what the
    # nMix=3 case below reproduces.
    threePopSplit <- function() {
      ini({
        tka   <- log(1.5)
        tcl1  <- log(1.0)
        tcl2  <- log(3.0)
        tcl3  <- log(6.0)
        tv    <- log(20)
        p1    <- 0.3
        p2    <- 0.4
        eta.cl1 ~ 0.01
        eta.cl2 ~ 0.01
        eta.cl3 ~ 0.01
        eta.v  ~ 0.01
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka)
        cl <- mix(exp(tcl1 + eta.cl1), p1, exp(tcl2 + eta.cl2), p2, exp(tcl3 + eta.cl3))
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .ui3 <- rxode2::rxode2(threePopSplit)
    .mkSaem3 <- function(mixProb) {
      .fixef <- setNames(rep(0.1, length(.ui3$saemParamsToEstimate)), .ui3$saemParamsToEstimate)
      .obj <- list(Plambda = .fixef,
                   resMat = matrix(rep(0.05, 4), nrow = 1),
                   mixProb = mixProb)
      class(.obj) <- "saemFit"
      .obj
    }

    # Two components each estimated near 0.7: individually valid probabilities,
    # but their sum (1.4) is >= 1 -- inconsistent/non-identifiable, requires
    # rescaling, and should warn.
    envCollapsed <- new.env()
    envCollapsed$ui <- .ui3
    envCollapsed$saem <- .mkSaem3(c(0.7, 0.7))
    expect_warning(
      nlmixr2est:::.getSaemTheta(envCollapsed),
      "collaps|mixture probabilit"
    )
    expect_equal(unname(envCollapsed$fullTheta[c("p1", "p2")]),
                 rep(0.7 / (1.4 + 1e-6), 2), tolerance = 1e-8)

    # Well-identified components: should stay silent, and values pass through
    # unchanged.
    envOk <- new.env()
    envOk$ui <- .ui3
    envOk$saem <- .mkSaem3(c(0.3, 0.4))
    expect_silent(nlmixr2est:::.getSaemTheta(envOk))
    expect_equal(unname(envOk$fullTheta[c("p1", "p2")]), c(0.3, 0.4))

    # Single component very near 0 (e.g. 1e-8): the raw estimate itself is
    # collapsed against the boundary, so this must warn even though the
    # clamp delta (1e-6 - 1e-8) is tiny in absolute terms -- collapse is
    # detected against the raw estimate, not the clamp/raw difference.
    one.compartment.mix <- function() {
      ini({
        tka <- log(1.5)
        tcl1 <- log(1.0)
        tcl2 <- log(5.0)
        tv <- log(20)
        p1 <- 0.5
        eta.cl ~ 0.01
        eta.v ~ 0.01
        eta.ka ~ 0.01
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .ui2 <- rxode2::rxode2(one.compartment.mix)
    .mkSaem2 <- function(mixProb) {
      .fixef <- setNames(c(log(1.5), log(1.0), log(5.0), log(20), 0.1),
                         .ui2$saemParamsToEstimate)
      .obj <- list(Plambda = .fixef,
                   resMat = matrix(rep(0.05, 4), nrow = 1),
                   mixProb = mixProb)
      class(.obj) <- "saemFit"
      .obj
    }
    envTinyBoundary <- new.env()
    envTinyBoundary$ui <- .ui2
    envTinyBoundary$saem <- .mkSaem2(1e-8)
    expect_warning(
      nlmixr2est:::.getSaemTheta(envTinyBoundary),
      "collaps|mixture probabilit"
    )
    expect_equal(unname(envTinyBoundary$fullTheta["p1"]), 1e-6, tolerance = 1e-8)

    # Well-identified two-component case (0.5): should stay silent.
    envOk2 <- new.env()
    envOk2$ui <- .ui2
    envOk2$saem <- .mkSaem2(0.5)
    expect_silent(nlmixr2est:::.getSaemTheta(envOk2))
    expect_equal(unname(envOk2$fullTheta["p1"]), 0.5)
  })

  test_that("test SAEM mixture model estimation", {

    # rxWithSeed pins BOTH the R and rxode2 RNG for the data sim and restores
    # them afterward, so this test neither depends on nor leaks seed state.
    rxode2::rxWithSeed(42, {
    n_subj <- 30
    sub_pop <- rbinom(n_subj, 1, 0.6) + 1 # 1 or 2
    cl_sim <- ifelse(sub_pop == 1, 1.2, 6.0)

    sim_data <- do.call(rbind, lapply(1:n_subj, function(i) {
      subj_cl <- cl_sim[i]
      times <- c(0.5, 1, 2, 4, 8, 12, 24)
      ka_val <- 1.5
      v_val <- 24.0
      k_val <- subj_cl / v_val
      cp <- 100 * ka_val / (v_val * (ka_val - k_val)) * (exp(-k_val * times) - exp(-ka_val * times)) + rnorm(length(times), 0, 0.05)
      cp[cp < 0] <- 0
      data.frame(
        ID = i,
        TIME = c(0, times),
        AMT = c(100, rep(0, length(times))),
        EVID = c(1, rep(0, length(times))),
        DV = c(0, cp),
        CMT = c(1, rep(2, length(times)))
      )
    }))
    })

    one.compartment.mix <- function() {
      ini({
        tka <- log(1.5)
        tcl1 <- log(1.0)
        tcl2 <- log(5.0)
        tv <- log(20)
        p1 <- 0.5
        eta.cl ~ 0.01
        eta.v ~ 0.01
        eta.ka ~ 0.01
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }

    # Run SAEM estimation with mixture
    fit_saem <- expect_error(
      .nlmixr(one.compartment.mix, sim_data, est="saem",
              saemControl(print = 0, seed = 1234, nBurn = 5, nEm = 5,
                          calcTables = FALSE, covMethod = 0L)),
      NA
    )

    # Check that estimated mixture parameter p1 exists
    expect_true("p1" %in% names(fit_saem$theta))

    # Check that estimated parameters are accessible
    expect_true(is.numeric(fit_saem$theta["p1"]))

    # Test SAEM mixture model with split ETAs and covariance calculation
    twoPopSplit <- function() {
      ini({
        tka   <- log(1.5)
        tcl1  <- log(1.0)
        tcl2  <- log(5.0)
        tv    <- log(20)
        p1    <- 0.5
        eta.cl1 ~ 0.01
        eta.cl2 ~ 0.01
        eta.v  ~ 0.01
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka)
        cl <- mix(exp(tcl1 + eta.cl1), p1, exp(tcl2 + eta.cl2))
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    fit_saem_split <- expect_error(
      .nlmixr(twoPopSplit, sim_data, est="saem",
              saemControl(print = 0, seed = 1234, nBurn = 5, nEm = 5,
                          calcTables = FALSE, covMethod = "linFim")),
      NA
    )
    expect_true("p1" %in% names(fit_saem_split$theta))
    expect_true(!is.null(fit_saem_split$cov))
    expect_true("eta.cl" %in% rownames(fit_saem_split$omega))
    expect_true(!"eta.cl1" %in% rownames(fit_saem_split$omega))
    expect_true("eta.cl" %in% names(fit_saem_split$ranef))
    expect_true(!"eta.cl1" %in% names(fit_saem_split$ranef))

    # Test SAEM mixture model with 3 split ETAs and covariance calculation
    threePopSplit <- function() {
      ini({
        tka   <- log(1.5)
        tcl1  <- log(1.0)
        tcl2  <- log(3.0)
        tcl3  <- log(6.0)
        tv    <- log(20)
        p1    <- 0.3
        p2    <- 0.4
        eta.cl1 ~ 0.01
        eta.cl2 ~ 0.01
        eta.cl3 ~ 0.01
        eta.v  ~ 0.01
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka)
        cl <- mix(exp(tcl1 + eta.cl1), p1, exp(tcl2 + eta.cl2), p2, exp(tcl3 + eta.cl3))
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    fit_saem_split3 <- expect_error(
      .nlmixr(threePopSplit, sim_data, est="saem",
              saemControl(print = 0, seed = 1234, nBurn = 5, nEm = 5,
                          calcTables = FALSE, covMethod = "linFim")),
      NA
    )
    expect_true("p1" %in% names(fit_saem_split3$theta))
    expect_true("p2" %in% names(fit_saem_split3$theta))
    expect_true(!is.null(fit_saem_split3$cov))
    expect_true("eta.cl" %in% rownames(fit_saem_split3$omega))
    expect_true(!"eta.cl1" %in% rownames(fit_saem_split3$omega))
    expect_true("eta.cl" %in% names(fit_saem_split3$ranef))
    expect_true(!"eta.cl1" %in% names(fit_saem_split3$ranef))

    # Test SAEM mixture model with split ETAs and bounded parameters
    # to verify that back-transformations are correctly applied (not NaN)
    twoPopSplitBounded <- function() {
      ini({
        tka   <- log(1.5)
        tcl1  <- log(c(0.01, 1.0, 100))  # Bounded!
        tcl2  <- log(c(0.01, 5.0, 100))  # Bounded!
        tv    <- log(20)
        p1    <- 0.5
        eta.cl1 ~ 0.01
        eta.cl2 ~ 0.01
        eta.v  ~ 0.01
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka)
        cl <- mix(exp(tcl1 + eta.cl1), p1, exp(tcl2 + eta.cl2))
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    fit_saem_split_bounded <- expect_error(
      .nlmixr(twoPopSplitBounded, sim_data, est="saem",
              saemControl(print = 0, seed = 1234, nBurn = 5, nEm = 5,
                          calcTables = FALSE, covMethod = "linFim")),
      NA
    )
    # Test SAEM mixture model with nested expit/logit functions under mix()
    # (non-split ETAs, so covariance cannot be calculated; covMethod=0L)
    twoPopNestedExpit <- function() {
      ini({
        tka   <- log(1.2)
        tcl1  <- logit(1.0, 0.1, 200)
        tcl2  <- logit(5.0, 0.1, 200)
        tv    <- log(30)
        p1    <- 0.67
        eta.cl ~ 0.09
        eta.v  ~ 0.04
        add.sd <- 0.2
      })
      model({
        ka <- exp(tka)
        cl <- mix(expit(tcl1 + eta.cl, 0.1, 200), p1, expit(tcl2 + eta.cl, 0.1, 200))
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    fit_saem_nested <- expect_error(
      .nlmixr(twoPopNestedExpit, sim_data, est="saem",
              saemControl(print = 0, seed = 1234, nBurn = 5, nEm = 5,
                          calcTables = FALSE, covMethod = 0L)),
      NA
    )
    expect_true("p1" %in% names(fit_saem_nested$theta))

    # Test SAEM mixture model with nested expit/logit functions under mix() and split ETAs
    # (split ETAs, so covariance is supported; covMethod="linFim")
    twoPopNestedExpitSplit <- function() {
      ini({
        tka   <- log(1.2)
        tcl1  <- logit(1.0, 0.1, 200)
        tcl2  <- logit(5.0, 0.1, 200)
        tv    <- log(30)
        p1    <- 0.67
        eta.cl1 ~ 0.09
        eta.cl2 ~ 0.09
        eta.v  ~ 0.04
        add.sd <- 0.2
      })
      model({
        ka <- exp(tka)
        cl <- mix(expit(tcl1 + eta.cl1, 0.1, 200), p1, expit(tcl2 + eta.cl2, 0.1, 200))
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    fit_saem_nested_split <- expect_error(
      .nlmixr(twoPopNestedExpitSplit, sim_data, est="saem",
              saemControl(print = 0, seed = 1234, nBurn = 5, nEm = 5,
                          calcTables = FALSE, covMethod = "linFim")),
      NA
    )
    expect_true("p1" %in% names(fit_saem_nested_split$theta))
    expect_true(!is.null(fit_saem_nested_split$cov))
  })

  test_that("SAEM mixture components actually separate (not just 'no error')", {
    # Regression test for the mixProb collapse bug (all subjects landing on
    # one component regardless of seed).
    # rxWithSeed pins BOTH the R and rxode2 RNG for the data sim and restores
    # them afterward, so this test neither depends on nor leaks seed state.
    rxode2::rxWithSeed(42, {
    n_subj <- 30
    sub_pop <- rbinom(n_subj, 1, 0.6) + 1 # 1 or 2
    cl_sim <- ifelse(sub_pop == 1, 1.2, 6.0)

    sim_data <- do.call(rbind, lapply(1:n_subj, function(i) {
      subj_cl <- cl_sim[i]
      times <- c(0.5, 1, 2, 4, 8, 12, 24)
      ka_val <- 1.5
      v_val <- 24.0
      k_val <- subj_cl / v_val
      cp <- 100 * ka_val / (v_val * (ka_val - k_val)) * (exp(-k_val * times) - exp(-ka_val * times)) + rnorm(length(times), 0, 0.05)
      cp[cp < 0] <- 0
      data.frame(
        ID = i,
        TIME = c(0, times),
        AMT = c(100, rep(0, length(times))),
        EVID = c(1, rep(0, length(times))),
        DV = c(0, cp),
        CMT = c(1, rep(2, length(times)))
      )
    }))
    })

    one.compartment.mix <- function() {
      ini({
        tka <- log(1.5)
        tcl1 <- log(1.0)
        tcl2 <- log(5.0)
        tv <- log(20)
        p1 <- 0.5
        eta.cl ~ 0.01
        eta.v ~ 0.01
        eta.ka ~ 0.01
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }

    fit <- .nlmixr(one.compartment.mix, sim_data, est="saem",
                   saemControl(print = 0, seed = 1234, nBurn = 200, nEm = 300,
                               calcTables = TRUE, covMethod = 0L))

    # A full collapse puts every subject in one component regardless of seed.
    mixTab <- table(fit$mixNum$mixnum)
    expect_true(length(mixTab) > 1)
    expect_true(min(mixTab) >= 2)

    # mixnum should track the true subpopulation better than a coin flip.
    agree <- mean(fit$mixNum$mixnum[order(fit$mixNum$ID)] == sub_pop)
    expect_true(agree > 0.6 || (1 - agree) > 0.6) # allow for label swap
  })

  test_that("SAEM split-ETA fixed effects actually separate (not just 'no error')", {
    # Regression test for split-ETA (separate eta.cl1/eta.cl2) collapse: tcl1
    # and tcl2 landing on nearly the same value even with mixProb/responsibility
    # fine. Checks back-transformed estimates are near truth and separated.
    #
    # NB: label order (which true subpopulation is component 1 vs 2) is not
    # guaranteed and can vary by seed -- no way to bound/order tcl1 vs tcl2
    # without breaking mu-referencing. So this test checks the
    # *set* of recovered clearances against the true set, not tcl1
    # specifically against the EM truth.
    # rxWithSeed pins the data stream and restores the global RNG state afterwards
    # (no bare set.seed leak into the fit or downstream tests)
    rxode2::rxWithSeed(2024, {
    tclEm <- log(8); tclPm <- log(0.8); tv0 <- log(30); tka0 <- log(1.2)
    sigma <- 0.20; omegaCl <- 0.09; omegaV <- 0.04
    nEm <- 20; nPm <- 10; n <- nEm + nPm
    clTrue <- c(exp(tclEm + rnorm(nEm, 0, sqrt(omegaCl))),
                exp(tclPm + rnorm(nPm, 0, sqrt(omegaCl))))
    vTrue <- exp(tv0 + rnorm(n, 0, sqrt(omegaV)))
    kaTrue <- rep(exp(tka0), n)
    times <- c(0.5, 1, 2, 4, 6, 8, 12, 18, 24, 36, 48)
    modSim <- rxode2::rxode2({
      d/dt(depot)   <- -ka * depot
      d/dt(central) <- ka * depot - (cl / v) * central
      cp <- central / v
    })
    simRows <- vector("list", n)
    for (i in seq_len(n)) {
      ev <- rxode2::et(amt = 100, time = 0) |> rxode2::et(time = times)
      out <- rxode2::rxSolve(modSim, params = c(ka = kaTrue[i], cl = clTrue[i], v = vTrue[i]), events = ev)
      dv <- out$cp * exp(rnorm(length(times), 0, sigma))
      simRows[[i]] <- data.frame(ID = i, time = times, DV = pmax(dv, 1e-6), AMT = 0, EVID = 0)
    }
    doseRows <- data.frame(ID = seq_len(n), time = 0, DV = NA_real_, AMT = 100, EVID = 1)
    simData <- rbind(doseRows, do.call(rbind, simRows))
    simData <- simData[order(simData$ID, simData$time), ]
    })

    # Nonlinear-wrapped (bounded) mu-referencing: cl <- mix(expit(tcl+eta,...))
    twoPopBounded <- function() {
      ini({
        tka <- log(1.2)
        tcl1 <- logit(0.8, 0.1, 200)
        tcl2 <- logit(0.8, 0.1, 200)
        tv <- log(30)
        p1 <- 0.67
        eta.cl1 ~ 0.09
        eta.cl2 ~ 0.09
        eta.v ~ 0.04
        add.sd <- 0.2
      })
      model({
        ka <- exp(tka)
        cl <- mix(expit(tcl1 + eta.cl1, 0.1, 200), p1, expit(tcl2 + eta.cl2, 0.1, 200))
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    fitBounded <- .nlmixr(twoPopBounded, simData, est="saem",
                          saemControl(print = 0, seed = 99, nBurn = 200, nEm = 300,
                                      calcTables = FALSE, covMethod = 0L))
    clBounded <- sort(c(fitBounded$theta[["tcl1"]], fitBounded$theta[["tcl2"]]))
    clBounded <- rxode2::expit(clBounded, 0.1, 200)
    # Recovered clearances should be near the true 0.8/8 and clearly
    # separated, not both near one of the true values.
    expect_true(clBounded[1] < 2)
    expect_true(clBounded[2] > 4)

    # Truly mu-referenced (linear, unbounded): cl <- mix(tcl+eta, p1, tcl+eta)
    twoPopMuLinear <- function() {
      ini({
        tka <- log(1.2)
        tcl1 <- 8
        tcl2 <- 0.8
        tv <- log(30)
        p1 <- 0.67
        eta.cl1 ~ 0.09
        eta.cl2 ~ 0.09
        eta.v ~ 0.04
        add.sd <- 0.2
      })
      model({
        ka <- exp(tka)
        cl <- mix(tcl1 + eta.cl1, p1, tcl2 + eta.cl2)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    fitMuLinear <- .nlmixr(twoPopMuLinear, simData, est="saem",
                           saemControl(print = 0, seed = 99, nBurn = 200, nEm = 300,
                                       calcTables = FALSE, covMethod = 0L))
    clMuLinear <- sort(c(fitMuLinear$theta[["tcl1"]], fitMuLinear$theta[["tcl2"]]))
    expect_true(clMuLinear[1] < 2)
    expect_true(clMuLinear[2] > 4)
  })
})


