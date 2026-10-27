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
    .ctl <- function(...) foceiControl(print = 0, maxOuterIterations = 0L, calcTables = FALSE, ...)
    # a step budget enough at the initial estimates (no outer iterations) but too
    # few for some of the covariance step's perturbed solves; where it falls varies
    # with the rxode2 build and platform
    .rs <- NULL
    for (.m in seq(85L, 130L, by = 5L)) {
      .fit <- .nlmixr(
        one.compartment,
        nlmixr2data::theo_sd,
        "focei",
        .ctl(rxControl = rxode2::rxControl(maxsteps = .m))
      )
      if (any(grepl("during the optimization", .fit$runInfo, fixed = TRUE))) {
        next
      }
      if (any(grepl("covariance step: atol/rtol", .fit$runInfo, fixed = TRUE))) {
        .rs <- .fit
        break
      }
    }
    if (is.null(.rs)) {
      skip("no maxsteps loosens only the covariance step's solves")
    }
    .none <- .nlmixr(
      one.compartment,
      nlmixr2data::theo_sd,
      "focei",
      .ctl(rxControl = rxode2::rxControl(maxsteps = .m), covMethod = "")
    )
    expect_false(any(grepl("tolerances", .none$runInfo, fixed = TRUE)))
    # the fit's own (and its tables') tolFactor; the covariance step's loosening
    # was reported as every subject's (316)
    expect_equal(.rs$env$tolFactor, .none$env$tolFactor)
  })

  test_that("a flat parameter needs no zero-gradient fallback; covSolveTol keeps tolFactor", {
    # A covariate coefficient whose covariate is 0 everywhere: the objective is
    # flat in tz.
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
    # with the ETAs held, the pooled gradient the S matrix falls back on is
    # exactly 0 in tz; that is a value there, not a zero the optimizer had replaced
    .held <- .nlmixr(flat, d, "focei", .ctl(maxInnerIterations = 0L))
    expect_false(any(grepl("zero gradient", .held$runInfo, fixed = TRUE)))

    # covSolveTol resets every subject's factor to 1 for the covariance solves;
    # a subject the estimation loosened kept 1 afterwards.  The step budget too
    # few for a subject at the fit's tolerance varies with the rxode2 build
    .none80 <- NULL
    for (.m in c(80L, 60L, 40L, 30L, 20L)) {
      .ctl80 <- function(...) {
        .ctl(rxControl = rxode2::rxControl(maxsteps = .m), ...)
      }
      .fit <- .nlmixr(one.compartment, nlmixr2data::theo_sd, "focei", .ctl80(covMethod = ""))
      if (any(.fit$env$tolFactor > 1)) {
        .none80 <- .fit
        break
      }
    }
    if (is.null(.none80)) {
      skip("no maxsteps loosens a subject during the estimation")
    }
    .tol80 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, "focei", .ctl80(covSolveTol = 1e-6))
    expect_equal(.tol80$env$tolFactor, .none80$env$tolFactor)
  })
})
