# Fit-based checks for the shared solve-pool registry (src/odeSwap.cpp).
#
# Split out of test-odeswap.R: these compile several models each, and enough
# model churn in the push/PR subset evicts entries from rxode2's model cache,
# which then breaks code generation in a LATER test file (seen as "user function
# 'expit' failed to produce code that could be parsed" in test-nlm.R).
#
# Weekly-batched via .slowBatches in tests/testthat.R -- do NOT add skip_on_ci().

nmTest({
  test_that("a real impmap model reaches the widest-lhs-is-not-the-pool case", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    # est="impmap" gives every residual-error theta a d(V)/d(theta) lhs column but
    # NO sensitivity state, so with more residual parameters than etas the
    # theta-sensitivity model is NARROWER in neq yet WIDER in lhs than the inner
    # model.  That is the configuration the private lhs buffer exists for; if this
    # ever stops holding, the scratch path is untested rather than unnecessary.
    m <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- fix(3.45)
        eta.ka ~ 0.6; eta.cl ~ 0.3
        add.sd <- 0.7; prop.sd <- 0.1; lambda <- 1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)
      })
    }
    ## The odeSwap counters are PROCESS-cumulative and nothing clears them, so a
    ## fit earlier in this worker leaves its arming behind.  Assert the DELTA
    ## across THIS fit; the absolute value only reads as 0 when the file happens
    ## to run before anything else that arms an override, which made it pass or
    ## fail on test-file scheduling rather than on behavior.
    .b <- .odeSwapInfo()
    .fit0 <- suppressWarnings(suppressMessages(
      nlmixr2(
        m,
        nlmixr2data::theo_sd,
        "impmap",
        impmapControl(
          print = 0L,
          nIter = 1L,
          isample = 50L,
          calcTables = FALSE,
          # combSens=TRUE (the default since #958) carries the theta
          # sensitivities on the INNER model, so no dedicated
          # theta-sensitivity peer is built at all.  This test is
          # about that peer, so ask for the two-model path.
          combSens = FALSE
        )
      )
    ))
    # Read the fit's OWN captured layout, not the process-global registry: the
    # registry describes the most recent registration, and impmap's post-fit
    # objective recompute runs a nested focei fit that re-registers the slots.
    i <- .fit0$env$odeSwapInfo
    ts <- i$models[i$models$name %in% "thetaSens", ]
    inr <- i$models[i$models$name %in% "inner", ]
    expect_identical(nrow(ts), 1L)
    expect_lt(ts$neq, inr$neq) # fewer states ...
    expect_gt(ts$nlhs, inr$nlhs) # ... but a wider lhs
    expect_identical(i$poolName, "inner")
    expect_false(ts$sizesPool)
    expect_true(i$needsScratch)
    expect_identical(i$scratchNlhs, ts$nlhs)
    # and the private buffer was actually taken during the fit -- without this the
    # test would still pass if OdeSwapScope silently handed back rxode2's slice
    expect_gt(i$scratchUsedN - .b$scratchUsedN, 0)
    # Inner solves now arm their own state-stride scope, so the process-wide
    # override counter cannot isolate the thetaSens solve.  scratchUsedN above
    # proves the private lhs buffer this test is about was actually used.
    expect_identical(i$scratchResizeN - .b$scratchResizeN, 0) # the plan sized it correctly up front
  })

  test_that("a fit does not inherit the previous fit's registered peers", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    # The registry is a global that outlives a fit.  A later fit with no
    # theta-sensitivity model of its own must NOT size its pool from the previous
    # fit's -- that is what resetting _impPoolModel around foceiSetup_ used to
    # prevent, and it shows up as "focei$rxInv needs to be of class
    # 'rxSymInvCholEnv'" on the following fit rather than as anything local.
    one <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    d <- nlmixr2data::theo_sd
    ctl <- foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L, calcTables = FALSE)
    ref <- suppressWarnings(suppressMessages(nlmixr2(one, d, "focei", ctl)))
    refPool <- .odeSwapInfo()
    expect_identical(refPool$poolName, "inner")

    # an impmap fit registers a theta-sensitivity peer that sizes the pool
    inv <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- fix(3.45)
        eta.ka ~ 0.6; eta.cl ~ 0.3
        add.sd <- 0.7; prop.sd <- 0.1; lambda <- 1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd) + prop(prop.sd) + boxCox(lambda)
      })
    }
    .fitInv <- suppressWarnings(suppressMessages(
      nlmixr2(
        inv,
        d,
        "impmap",
        impmapControl(
          print = 0L,
          nIter = 1L,
          isample = 50L,
          calcTables = FALSE,
          # combSens=TRUE (the default since #958) carries the theta
          # sensitivities on the INNER model, so no dedicated
          # theta-sensitivity peer is built at all.  This test is
          # about that peer, so ask for the two-model path.
          combSens = FALSE
        )
      )
    ))
    # the fit's own captured layout -- see the note above
    expect_true(.fitInv$env$odeSwapInfo$models$loaded[3]) # thetaSens registered

    # ... and the next plain fit must be unaffected by it
    after <- suppressWarnings(suppressMessages(nlmixr2(one, d, "focei", ctl)))
    expect_identical(.odeSwapInfo()$poolName, "inner")
    expect_identical(after$objf, ref$objf)
    expect_identical(unname(as.numeric(after$theta)), unname(as.numeric(ref$theta)))
  })

  test_that("a pinned inner override does not leak into the next fit", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    # est="impmap" pins every subject's effective state count for the whole fit
    # (the pool is sized for the larger theta-sensitivity model).  That pin used
    # to be released only by the impmap driver, so the fast-ll and vae paths left
    # it set on the shared solve structure for whatever ran next.  The release now
    # lives in rxOptionsFreeFocei(), which runs at BOTH setup start and teardown.
    d <- nlmixr2data::theo_sd
    # 1-compartment reference
    one <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
            eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              d/dt(depot) <- -ka * depot
              d/dt(center) <- ka * depot - cl / v * center
              cp <- center / v; cp ~ add(add.sd) })
    }
    # the pinning fit, with MORE states -- a stale pin only bites when the next
    # fit's op->neq differs from the one that was pinned
    pin <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- fix(3.45); tq <- 0.1; tvp <- 1
            eta.ka ~ 0.6; eta.cl ~ 0.3
            add.sd <- 0.7; prop.sd <- 0.1; lambda <- 1 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
              q <- exp(tq); vp <- exp(tvp)
              d/dt(depot) <- -ka * depot
              d/dt(center) <- ka * depot - cl / v * center - q/v * center + q/vp * peri
              d/dt(peri) <- q/v * center - q/vp * peri
              cp <- center / v
              cp ~ add(add.sd) + prop(prop.sd) + boxCox(lambda) })
    }
    ctl <- foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L, calcTables = FALSE)
    ref <- suppressWarnings(suppressMessages(nlmixr2(one, d, "focei", ctl)))
    expect_false(.odeSwapInfo()$pinned)

    suppressWarnings(suppressMessages(
      nlmixr2(pin, d, "impmap", impmapControl(print = 0L, nIter = 1L, isample = 50L, calcTables = FALSE))
    ))
    # the pin is released by the time the fit returns
    expect_false(.odeSwapInfo()$pinned)
    expect_true(all(.odeSwapInfo()$activeOverride == -1L))

    after <- suppressWarnings(suppressMessages(nlmixr2(one, d, "focei", ctl)))
    expect_identical(after$objf, ref$objf)
    expect_identical(unname(as.numeric(after$theta)), unname(as.numeric(ref$theta)))
  })

  test_that("a pooled multi-endpoint model re-bases CMT and stays correct", {
    skip_on_cran()
    # rxode2 normalizes CMT inside each compiled model with THAT model's own
    # sensitivity-compartment count (`#define _CMT ... CMT - nSens`), so one translated
    # event table cannot serve peers of different sensitivity depth.  Pooled, the inner
    # model was handed the augmented model's basis and computed 63 - 2 = 61, matching no
    # endpoint: rx_pred_, rx_r_, d(f)/d(eta) and rx_yj_ all evaluated to 0, the EBEs
    # collapsed to ~0, and yj = 0 silently log-transformed DV.  OdeSwapCmtScope re-bases
    # the CMT covariate per solving model; this pins that it holds.
    #
    # Asserts the MECHANISM (which model sized the pool) as well as the result: a
    # matching objective alone would also pass if fast=TRUE had quietly stopped
    # being fast for an unrelated reason.
    me <- function() {
      ini({ tka <- 0.5; tcl <- -2; tv <- 2; tbase <- 1; tkout <- -1
            add.pk <- 1; add.pd <- 0.5; eta.cl ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
              d/dt(depot) <- -ka * depot
              d/dt(center) <- ka * depot - cl / v * center
              cp <- center / v
              pd <- exp(tbase) * cp / (exp(tkout) + cp)
              cp ~ add(add.pk) | cp
              pd ~ add(add.pd) | pd })
    }
    withr::local_seed(7)
    mk <- function(i, mu) {
      o <- rbind(
        data.frame(ID = i, TIME = c(.5, 1, 2, 4, 8), EVID = 0, AMT = 0, DV = abs(stats::rnorm(5, mu, 1)), CMT = "cp"),
        data.frame(
          ID = i,
          TIME = c(.75, 1.5, 3, 6, 10),
          EVID = 0,
          AMT = 0,
          DV = abs(stats::rnorm(5, 2, .3)),
          CMT = "pd"
        )
      )
      rbind(data.frame(ID = i, TIME = 0, EVID = 101, AMT = 100, DV = NA, CMT = "cp"), o[order(o$TIME), ])
    }
    d <- rbind(mk(1, 8), mk(2, 3))
    # one thread: the parallel inner optimizer is not bitwise reproducible, and this
    # compares EBEs, not just the objective
    .th <- rxode2::rxCores()
    on.exit(rxode2::setRxThreads(.th), add = TRUE)
    rxode2::setRxThreads(1L)
    ctl <- function(fast) {
      foceiControl(
        print = 0L,
        covMethod = "",
        sigdig = 4,
        fast = fast,
        maxOuterIterations = 0L,
        maxInnerIterations = 200L,
        calcTables = FALSE
      )
    }
    ref <- suppressWarnings(suppressMessages(nlmixr2(me, d, "focei", ctl(FALSE))))
    fast <- suppressWarnings(suppressMessages(nlmixr2(me, d, "focei", ctl(TRUE))))
    # the augmented model DOES size the pool here -- the point is that the CMT
    # re-base keeps the objective right anyway
    expect_identical(.odeSwapInfo()$poolName, "outer")
    ## not expect_identical: fast=TRUE still evaluates the posthoc gradient, so the
    ## last digits move.  The bug this pins was worth ~900 objective units, not 1e-12.
    expect_equal(as.numeric(fast$objf), as.numeric(ref$objf), tolerance = 1e-8)
    expect_equal(as.numeric(fast$eta[[2]]), as.numeric(ref$eta[[2]]), tolerance = 1e-6)
    # and the ETAs are genuinely conditional, not the collapsed ~0 the shift produced
    expect_gt(max(abs(as.numeric(ref$eta[[2]]))), 1e-3)
  })

  test_that("an infusion model keeps its inner state stride in the outer pool", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    wbc <- function() {
      ini({
        log_CIRC0 <- log(7.21); log_MTT <- log(124)
        log_SLOPU <- log(28.9); log_GAMMA <- log(0.239)
        eta.CIRC0 ~ 0.1; eta.MTT ~ 0.03; eta.SLOPU ~ 0.2
        prop.err <- 10
      })
      model({
        CIRC0 <- exp(log_CIRC0 + eta.CIRC0)
        MTT <- exp(log_MTT + eta.MTT)
        SLOPU <- exp(log_SLOPU + eta.SLOPU)
        GAMMA <- exp(log_GAMMA)
        CL <- CLI; V1 <- V1I; V2 <- V2I; Q <- 204
        CONC <- centr / V1
        KTR <- 4 / MTT
        EDRUG <- 1 - SLOPU * CONC
        FDBK <- (CIRC0 / circ)^GAMMA
        CIRC <- circ
        prol(0) <- CIRC0; tr1(0) <- CIRC0; tr2(0) <- CIRC0
        tr3(0) <- CIRC0; circ(0) <- CIRC0
        d/dt(centr) <- periph * Q / V2 - centr * (CL / V1 + Q / V1)
        d/dt(periph) <- centr * Q / V1 - periph * Q / V2
        d/dt(prol) <- KTR * prol * EDRUG * FDBK - KTR * prol
        d/dt(tr1) <- KTR * prol - KTR * tr1
        d/dt(tr2) <- KTR * tr1 - KTR * tr2
        d/dt(tr3) <- KTR * tr2 - KTR * tr3
        d/dt(circ) <- KTR * tr3 - KTR * circ
        CIRC ~ prop(prop.err)
      })
    }
    d <- nlmixr2data::wbcSim[nlmixr2data::wbcSim$ID == 1, ]
    ctl <- function(fast) {
      foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L,
                   calcTables = FALSE, fast = fast)
    }
    ref <- suppressWarnings(suppressMessages(nlmixr2(wbc, d, "focei", ctl(FALSE))))
    fast <- suppressWarnings(suppressMessages(nlmixr2(wbc, d, "focei", ctl(TRUE))))
    expect_identical(.odeSwapInfo()$poolName, "outer")
    expect_equal(as.numeric(fast$objf), as.numeric(ref$objf), tolerance = 1e-7)
    expect_equal(as.numeric(fast$eta[[2]]), as.numeric(ref$eta[[2]]), tolerance = 1e-6)
  })

  test_that("the nonlinear nimo model keeps its inner solve and analytic sensitivities", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    nimo <- function() {
      ini({
        tcl <- log(0.001); tv1 <- log(1.45); tQ <- log(0.004); tv2 <- log(44)
        tkss <- log(12); tkint <- log(0.3); tksyn <- log(1); tkdeg <- log(7)
        eta.cl ~ 2; eta.v1 ~ 2; eta.kss ~ 2; add.err <- 10
      })
      model({
        cl <- exp(tcl + eta.cl); v1 <- exp(tv1 + eta.v1)
        Q <- exp(tQ); v2 <- exp(tv2); kss <- exp(tkss + eta.kss)
        kint <- exp(tkint); ksyn <- exp(tksyn); kdeg <- exp(tkdeg)
        k <- cl / v1; k12 <- Q / v1; k21 <- Q / v2
        eff(0) <- ksyn / kdeg
        conc <- 0.5 * (central / v1 - eff - kss) +
          0.5 * sqrt((central / v1 - eff - kss)^2 + 4 * kss * central / v1)
        d/dt(central) <- -(k + k12) * conc * v1 + k21 * peripheral -
          kint * eff * conc * v1 / (kss + conc)
        d/dt(peripheral) <- k12 * conc * v1 - k21 * peripheral
        d/dt(eff) <- ksyn - kdeg * eff - (kint - kdeg) * conc * eff / (kss + conc)
        IPRED <- log(conc)
        IPRED ~ add(add.err)
      })
    }
    ctl <- function(fast) {
      foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L,
                   calcTables = FALSE, fast = fast)
    }
    ref <- suppressWarnings(suppressMessages(
      nlmixr2(nimo, nlmixr2data::nimoData, "focei", ctl(FALSE))))
    fast <- suppressWarnings(suppressMessages(
      nlmixr2(nimo, nlmixr2data::nimoData, "focei", ctl(TRUE))))
    expect_identical(.odeSwapInfo()$poolName, "outer")
    expect_equal(as.numeric(fast$objf), as.numeric(ref$objf), tolerance = 1e-6)
    expect_equal(as.matrix(fast$eta[-1]), as.matrix(ref$eta[-1]), tolerance = 1e-5)
  })

  test_that("a parameter-dependent initial condition converges in a pooled fast fit", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    m <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.cl ~ 0.3; add.sd <- 0.7 })
      model({
        ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        center(0) <- v * 0.001
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    ctl <- function(fast) {
      foceiControl(print = 0L, covMethod = "", calcTables = FALSE, fast = fast)
    }
    ref <- suppressWarnings(suppressMessages(
      nlmixr2(m, nlmixr2data::theo_sd, "focei", ctl(FALSE))))
    fast <- suppressWarnings(suppressMessages(
      nlmixr2(m, nlmixr2data::theo_sd, "focei", ctl(TRUE))))
    expect_equal(as.numeric(fast$objf), as.numeric(ref$objf), tolerance = 1e-4)
  })
})
