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
    suppressWarnings(suppressMessages(
      nlmixr2(m, nlmixr2data::theo_sd, "impmap",
              impmapControl(print = 0L, nIter = 1L, isample = 50L, calcTables = FALSE))))
    i <- .odeSwapInfo()
    ts <- i$models[i$models$name %in% "thetaSens", ]
    inr <- i$models[i$models$name %in% "inner", ]
    expect_identical(nrow(ts), 1L)
    expect_lt(ts$neq, inr$neq)     # fewer states ...
    expect_gt(ts$nlhs, inr$nlhs)   # ... but a wider lhs
    expect_identical(i$poolName, "inner")
    expect_false(ts$sizesPool)
    expect_true(i$needsScratch)
    expect_identical(i$scratchNlhs, ts$nlhs)
    # and the private buffer was actually taken during the fit -- without this the
    # test would still pass if OdeSwapScope silently handed back rxode2's slice
    expect_gt(i$overrideArmedN, 0)
    expect_gt(i$scratchUsedN, 0)
    expect_identical(i$scratchResizeN, 0)   # the plan sized it correctly up front
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
    ctl <- foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L,
                        calcTables = FALSE)
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
    suppressWarnings(suppressMessages(
      nlmixr2(inv, d, "impmap",
              impmapControl(print = 0L, nIter = 1L, isample = 50L, calcTables = FALSE))))
    expect_true(odeSwapInfo_()$models$loaded[3])   # thetaSens registered

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
    ctl <- foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L,
                        calcTables = FALSE)
    ref <- suppressWarnings(suppressMessages(nlmixr2(one, d, "focei", ctl)))
    expect_false(.odeSwapInfo()$pinned)

    suppressWarnings(suppressMessages(
      nlmixr2(pin, d, "impmap",
              impmapControl(print = 0L, nIter = 1L, isample = 50L, calcTables = FALSE))))
    # the pin is released by the time the fit returns
    expect_false(.odeSwapInfo()$pinned)
    expect_true(all(.odeSwapInfo()$activeOverride == -1L))

    after <- suppressWarnings(suppressMessages(nlmixr2(one, d, "focei", ctl)))
    expect_identical(after$objf, ref$objf)
    expect_identical(unname(as.numeric(after$theta)),
                     unname(as.numeric(ref$theta)))
  })
})
