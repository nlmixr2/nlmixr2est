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

  # Shared base fit (with default calcTables so the result is nlmixr2FitData).
  # NPDE is not requested here so addNpde tests can add it themselves.
  fit <- .nlmixr(
    one.compartment,
    nlmixr2data::theo_sd,
    est = "focei",
    control = foceiControl(print = 0, covMethod = "")
  )

  test_that("tolFactor output: fit$env$tolFactor is a named numeric vector with one entry per subject >= 1", {
    tf <- fit$env$tolFactor
    nsubj <- length(unique(nlmixr2data::theo_sd$ID))
    expect_true(is.numeric(tf))
    expect_equal(length(tf), nsubj)
    expect_true(all(tf >= 1))
  })

  test_that("tolFactor output is used in addCwres: CWRES are finite", {
    # addCwres returns early if CWRES already exists (computed during default
    # table calculation), so the check still exercises the tolFactor code path.
    fitC <- suppressMessages(addCwres(fit))
    expect_true("CWRES" %in% names(fitC))
    expect_true(all(is.finite(fitC$CWRES)))
  })

  test_that("tolFactor output is used in addNpde: NPDE are finite", {
    fitN <- suppressMessages(
      addNpde(fit, table = tableControl(nsim = 50, seed = 42))
    )
    expect_true("NPDE" %in% names(fitN))
    expect_true(all(is.finite(fitN$NPDE)))
  })

  test_that("tolFactor output is used when npde is requested at fit time via tableControl", {
    fitNpde <- .nlmixr(
      one.compartment,
      nlmixr2data::theo_sd,
      est = "focei",
      control = foceiControl(print = 0, covMethod = ""),
      table = tableControl(npde = TRUE, nsim = 50, seed = 42)
    )
    expect_true("NPDE" %in% names(fitNpde))
    expect_true(all(is.finite(fitNpde$NPDE)))
    tf <- fitNpde$env$tolFactor
    nsubj <- length(unique(nlmixr2data::theo_sd$ID))
    expect_equal(length(tf), nsubj)
    expect_true(all(tf >= 1))
  })

  test_that("tolFactor input: rxControl(tolFactor=5) is stored in the fit and used by the ODE solver", {
    fit5 <- .nlmixr(
      one.compartment,
      nlmixr2data::theo_sd,
      est = "focei",
      control = foceiControl(print = 0, covMethod = "", rxControl = rxode2::rxControl(tolFactor = 5))
    )
    expect_s3_class(fit5, "nlmixr2FitData")

    # The input tolFactor is preserved in the stored rxControl and therefore
    # passed to every rxSolve_ call made with these options.
    expect_equal(fit5$foceiControl$rxControl$tolFactor, 5)

    # The per-subject output tolFactor accumulated during FOCEI optimization
    # is a valid numeric vector (>= 1); it is separate from the rxControl
    # tolFactor, which seeds the ODE-solver atol/rtol for the initial setup solve.
    tf <- fit5$env$tolFactor
    nsubj <- length(unique(nlmixr2data::theo_sd$ID))
    expect_equal(length(tf), nsubj)
    expect_true(all(tf >= 1))

    # addCwres and addNpde use fit$env$tolFactor for post-fit solves; verify
    # they work even when the fit was produced with a non-default rxControl tolFactor.
    fit5C <- suppressMessages(addCwres(fit5))
    expect_true(all(is.finite(fit5C$CWRES)))

    fit5N <- suppressMessages(
      addNpde(fit5, table = tableControl(nsim = 50, seed = 42))
    )
    expect_true(all(is.finite(fit5N$NPDE)))
  })

  test_that("the covariance step leaves the fit's tolFactor and its warnings alone", {
    # A covariate coefficient whose covariate is 0 everywhere: the objective is
    # flat in tz, so the covariance step's search for its step size probes far
    # enough to need looser tolerances.  The estimation (no outer iterations)
    # never does.
    flat <- function() {
      ini({
        tka <- log(1.5)
        tcl <- log(2.7)
        tv <- log(31.5)
        tz <- 0.1
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + tz * Z)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    d <- nlmixr2data::theo_sd
    d$Z <- 0
    .ctl <- function(...) foceiControl(print = 0, maxOuterIterations = 0L, calcTables = FALSE, ...)
    .none <- .nlmixr(flat, d, "focei", .ctl(covMethod = ""))
    .rs <- .nlmixr(flat, d, "focei", .ctl())
    expect_false(any(grepl("tolerances", .none$runInfo, fixed = TRUE)))
    # the fit's own (and its tables') tolFactor; the covariance step's loosening
    # was reported as every subject's (316)
    expect_equal(.rs$env$tolFactor, .none$env$tolFactor)
    # and the loosening is the covariance step's, not the optimization's
    expect_true(any(grepl("covariance step: atol/rtol increased", .rs$runInfo, fixed = TRUE)))
    expect_false(any(grepl("during the optimization", .rs$runInfo, fixed = TRUE)))
    # with the ETAs held, the pooled gradient the S matrix falls back on is
    # exactly 0 in tz; that is a value there, not a zero the optimizer had replaced
    .held <- .nlmixr(flat, d, "focei", .ctl(maxInnerIterations = 0L))
    expect_false(any(grepl("zero gradient", .held$runInfo, fixed = TRUE)))

    # covSolveTol resets every subject's factor to 1 for the covariance solves;
    # a subject the estimation loosened (maxsteps = 80 is too few for one of them
    # at the fit's tolerance) kept 1 afterwards
    .ctl80 <- function(...) {
      .ctl(rxControl = rxode2::rxControl(maxsteps = 80L), ...)
    }
    .none80 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, "focei", .ctl80(covMethod = ""))
    expect_true(any(.none80$env$tolFactor > 1))
    .tol80 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, "focei", .ctl80(covSolveTol = 1e-6))
    expect_equal(.tol80$env$tolFactor, .none80$env$tolFactor)
  })
})
