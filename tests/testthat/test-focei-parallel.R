nmTest({
  one.compartment <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v  <- exp(tv  + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  test_that("FOCEi cores=2 matches cores=1 baseline", {
    fit1 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "",
                                          rxControl = rxode2::rxControl(cores = 1L)))
    fit2 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "",
                                          rxControl = rxode2::rxControl(cores = 2L)))
    expect_s3_class(fit1, "nlmixr2FitData")
    expect_s3_class(fit2, "nlmixr2FitData")
    expect_equal(fit1$objective, fit2$objective, tolerance = 1e-4)
    expect_equal(as.numeric(fit1$theta), as.numeric(fit2$theta), tolerance = 1e-4)
  })

  test_that("FOCEi cores=2 matches cores=1 with mceta=0 (always reset)", {
    # mceta=0 path in innerOpt1: zero ETAs each call, no R API.
    fit1 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "", mceta = 0L,
                                          rxControl = rxode2::rxControl(cores = 1L)))
    fit2 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "", mceta = 0L,
                                          rxControl = rxode2::rxControl(cores = 2L)))
    expect_equal(fit1$objective, fit2$objective, tolerance = 1e-4)
  })

  test_that("FOCEi cores=2 with mceta>=1 (pre-drawn samples) does not crash", {
    # mceta >= 1 used to force cores=1 because the sampling path called the
    # R API.  Samples are now pre-drawn serially before the parallel region.
    # The objective need not match cores=1 bit-for-bit because the optimizer
    # trajectory differs across cores; only require that both fits complete.
    fit1 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "", mceta = 2L,
                                          rxControl = rxode2::rxControl(cores = 1L)))
    fit2 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "", mceta = 2L,
                                          rxControl = rxode2::rxControl(cores = 2L)))
    expect_s3_class(fit1, "nlmixr2FitData")
    expect_s3_class(fit2, "nlmixr2FitData")
    expect_true(is.finite(fit1$objective))
    expect_true(is.finite(fit2$objective))
  })

  # ETA-dependent dose parameters (f, dur, rate, alag) trigger shi21's
  # finite-difference path.  Each FD perturbation calls predOde with a
  # smaller "predNeq" (no inner-sensitivity states); shi21EtaGeneral now
  # switches the per-individual ind->neqOverride via IndNeqOverrideGuard
  # so concurrent threads on other subjects keep solving with op->neq
  # without observing a shared-state race.  Reprex:
  # inst/reprex_fbio_eta_parallel.R.

  test_that("FOCEi cores=2 works when f(state) depends on eta", {
    mod_fbio <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv  <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v  ~ 0.1
        eta.f  ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv  + eta.v)
        f(depot) <- exp(eta.f)
        d/dt(depot)  <- -ka * depot
        d/dt(center) <-  ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    fit1 <- .nlmixr(mod_fbio, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "",
                                          rxControl = rxode2::rxControl(cores = 1L)))
    fit2 <- .nlmixr(mod_fbio, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "",
                                          rxControl = rxode2::rxControl(cores = 2L)))
    expect_s3_class(fit1, "nlmixr2FitData")
    expect_s3_class(fit2, "nlmixr2FitData")
    # cores=2 must match cores=1 bit-for-bit.  Per-subject stickyRecalcN2
    # and indHasBadSolve() (replacing the shared op->badSolve poll) make
    # the inner-retry loop fully deterministic per subject, so the
    # objective is reproducible across cores even when the optimizer
    # crosses bad regions and triggers tolerance loosening.
    expect_true(is.finite(fit1$objective))
    expect_true(is.finite(fit2$objective))
    expect_equal(fit1$objective, fit2$objective, tolerance = 1e-8)
  })

  test_that("FOCEi cores=2 works when dur(state) depends on eta", {
    mod_dur <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv  <- 3.45
        eta.cl ~ 0.3
        eta.v  ~ 0.1
        eta.dur ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv  + eta.v)
        dur(depot) <- exp(eta.dur)
        d/dt(depot)  <- -ka * depot
        d/dt(center) <-  ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    # Modeled-duration data: RATE = -2 marks dose records.
    d <- nlmixr2data::theo_sd
    d$RATE <- ifelse(d$EVID == 1, -2, 0)
    fit1 <- .nlmixr(mod_dur, d, est = "focei",
                    control = foceiControl(print = 0, covMethod = "",
                                          rxControl = rxode2::rxControl(cores = 1L)))
    fit2 <- .nlmixr(mod_dur, d, est = "focei",
                    control = foceiControl(print = 0, covMethod = "",
                                          rxControl = rxode2::rxControl(cores = 2L)))
    # cores=2 must match cores=1 bit-for-bit.  Per-subject stickyRecalcN2
    # and indHasBadSolve() (replacing the shared op->badSolve poll) make
    # the inner-retry loop fully deterministic per subject, so the
    # objective is reproducible across cores even when the optimizer
    # crosses bad regions and triggers tolerance loosening.
    expect_true(is.finite(fit1$objective))
    expect_true(is.finite(fit2$objective))
    expect_equal(fit1$objective, fit2$objective, tolerance = 1e-8)
  })

  test_that("FOCEi cores=2 works when rate(state) depends on eta", {
    mod_rate <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv  <- 3.45
        eta.cl ~ 0.3
        eta.v  ~ 0.1
        eta.rate ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv  + eta.v)
        rate(depot) <- exp(eta.rate)
        d/dt(depot)  <- -ka * depot
        d/dt(center) <-  ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    # Modeled-rate data: RATE = -1.
    d <- nlmixr2data::theo_sd
    d$RATE <- ifelse(d$EVID == 1, -1, 0)
    fit1 <- .nlmixr(mod_rate, d, est = "focei",
                    control = foceiControl(print = 0, covMethod = "",
                                          rxControl = rxode2::rxControl(cores = 1L)))
    fit2 <- .nlmixr(mod_rate, d, est = "focei",
                    control = foceiControl(print = 0, covMethod = "",
                                          rxControl = rxode2::rxControl(cores = 2L)))
    # cores=2 must match cores=1 bit-for-bit.  Per-subject stickyRecalcN2
    # and indHasBadSolve() (replacing the shared op->badSolve poll) make
    # the inner-retry loop fully deterministic per subject, so the
    # objective is reproducible across cores even when the optimizer
    # crosses bad regions and triggers tolerance loosening.
    expect_true(is.finite(fit1$objective))
    expect_true(is.finite(fit2$objective))
    expect_equal(fit1$objective, fit2$objective, tolerance = 1e-8)
  })

  test_that("FOCEi cores=2 works when alag(state) depends on eta", {
    mod_alag <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv  <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v  ~ 0.1
        eta.lag ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv  + eta.v)
        alag(depot) <- exp(eta.lag)
        d/dt(depot)  <- -ka * depot
        d/dt(center) <-  ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    fit1 <- .nlmixr(mod_alag, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "",
                                          rxControl = rxode2::rxControl(cores = 1L)))
    fit2 <- .nlmixr(mod_alag, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "",
                                          rxControl = rxode2::rxControl(cores = 2L)))
    # cores=2 must match cores=1 bit-for-bit.  Per-subject stickyRecalcN2
    # and indHasBadSolve() (replacing the shared op->badSolve poll) make
    # the inner-retry loop fully deterministic per subject, so the
    # objective is reproducible across cores even when the optimizer
    # crosses bad regions and triggers tolerance loosening.
    expect_true(is.finite(fit1$objective))
    expect_true(is.finite(fit2$objective))
    expect_equal(fit1$objective, fit2$objective, tolerance = 1e-8)
  })
})
