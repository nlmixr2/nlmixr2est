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

    # Single component very near 0 (e.g. 1e-8): clamped to 1e-6, but the
    # change is tiny in absolute terms, so this does NOT cross the warning
    # threshold on its own -- documenting the boundary behavior explicitly.
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
    expect_silent(nlmixr2est:::.getSaemTheta(envTinyBoundary))
    expect_equal(unname(envTinyBoundary$fullTheta["p1"]), 1e-6, tolerance = 1e-8)

    # Well-identified two-component case (0.5): should stay silent.
    envOk2 <- new.env()
    envOk2$ui <- .ui2
    envOk2$saem <- .mkSaem2(0.5)
    expect_silent(nlmixr2est:::.getSaemTheta(envOk2))
    expect_equal(unname(envOk2$fullTheta["p1"]), 0.5)
  })

  test_that("test SAEM mixture model estimation", {

    set.seed(42)
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
})


