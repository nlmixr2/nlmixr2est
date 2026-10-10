nmTest({
  test_that("focei complex event info", {
    pheno <- function() {
      ini({
        tcl <- log(0.008) # typical value of clearance
        tv <-  log(0.6)   # typical value of volume
        max_dose <- 5
        ## var(eta.cl)
        eta.cl + eta.v ~ c(1,
                           0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
        # interindividual variability on clearance and volume
        add.err <- 0.1    # residual variability
      })
      model({
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        fest <- max_dose/DOSE
        if (fest > 1) fest <- 1 # error is here
        d/dt(A1) = - ke * A1
        f(A1) <- fest
        cp = A1 / v
        cp ~ add(add.err)
      })
    }

    f <- suppressMessages(pheno())
    expect_error(f$foceiModel, NA)
  })

  test_that("Inner test", {
    ev <- eventTable() |>
      add.sampling(c(
        95.99,
        119.99,
        143.99,
        144.25,
        144.5,
        144.75,
        145,
        145.5,
        146,
        146.5,
        147,
        148,
        150,
        152,
        156,
        160,
        164,
        167.99,
        191.99,
        215.99,
        216.25,
        216.5,
        216.75,
        217,
        217.5,
        218,
        218.5,
        219,
        220,
        222,
        224,
        228,
        232,
        236,
        240,
        252,
        264,
        276,
        288
      )) |>
      add.dosing(dose = 60000, start.time = 72, nbr.doses = 7, dosing.interval = 24)

    dv <- c(
      263.6,
      164.7,
      287.3,
      1248.7,
      1211.5,
      1017.7,
      1690.1,
      1029.8,
      890.7,
      598.4,
      1009.3,
      1159.8,
      742.2,
      724.6,
      728.2,
      509.7,
      243.1,
      259.9,
      242.2,
      281.4,
      1500.1,
      1281.4,
      1200.2,
      1378.8,
      1373.2,
      582.9,
      960.2,
      720.3,
      852.6,
      950.3,
      654.7,
      402.5,
      456,
      346.5,
      268.2,
      134.2,
      42.6,
      25.9,
      14.6
    )

    m1 <- function() {
      ini({
        tcl <- 1.6
        tv <- 4.5
        eta.cl ~ 0.1
        eta.v ~ 0.1
        prop.sd <- sqrt(0.1)
      })
      model({
        CL <- exp(tcl + eta.cl)
        V <- exp(tv + eta.v)
        C2 <- centr / V
        d/dt(centr) <- -CL * C2
        C2 ~ prop(prop.sd)
      })
    }

    w7 <- data.frame(ev$get.EventTable())
    w7$DV <- NA
    w7$DV[which(is.na(w7$amt))] <- dv
    w7$ID <- 1

    ETA <- matrix(c(-0.147736086922763, -0.294637022436797), ncol = 2)

    ## sigdig pinned: this asserts the CONVERGED inner objective to three
    ## decimals, and sigdig drives the solver tolerances (rtol = 10^-sigdig).
    ## At the sigdig=3 default (rtol=1e-3) the solve cannot resolve that third
    ## decimal -- it returns 418.9344 where sigdig 4/5/6 all give 418.9353 -- so
    ## the assertion failed on solver precision, not on the objective.  The
    ## reference value is unchanged; only the solve is now tight enough to
    ## reproduce it independently of whatever the default sigdig happens to be.
    fitPi <- .nlmixr(
      m1,
      w7,
      est = "focei",
      foceiControl(
        etaMat = ETA,
        maxOuterIterations = 0,
        maxInnerIterations = 0,
        covMethod = "",
        sigdig = 6
      )
    )

    expect_equal(418.935, round(fitPi$objective, 3))
  })

  test_that("boundary value is not triggered by bounds on both sides of zero (#318)", {
    one.compartment <- function() {
      ini({
        tka <- c(-6, -4, 2)
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)*100
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.sd)
      })
    }

    fit <- .nlmixr(one.compartment, theo_sd, est = "focei", control = list(print = 0))
    # SE being present indicates that the covariance matrix was estimated
    expect_true("SE" %in% names(fit$parFixedDf))

    # Also make sure that it correctly identifies mu-ref
    expect_equal(c("tka", "tcl", "tv", "add.sd"), row.names(fit$parFixedDf))
  })

  test_that("focei model with sine over a compound argument builds (#513)", {
    # A trig function whose argument is a compound expression divided by
    # something (here 2*3.14*(time-mtime1)/period) used to lose its argument in
    # the symengine round-trip, emitting sin()/cos() with no argument and
    # failing to compile ("too few arguments to function 'sin'").  Requires the
    # rxFromSE fix in rxode2; skip on an rxode2 that still drops the argument.
    skip_if(rxode2::rxFromSE("sin((a-b)/c)") == "sin()", "installed rxode2 predates the rxFromSE compound-argument fix")

    ehc <- function() {
      ini({
        tKa <- log(10)
        tCl <- log(93794.73)
        tV <- log(973551.9)
        tKemp <- log(30)
        tmtime1 <- log(1)
        tperiod <- 6
        add.err <- c(0, 0.01)
        prop.err <- c(0, 0.2)
        eta.Ka ~ 0.1
        eta.mtime1 ~ 0.1
      })
      model({
        Ka <- exp(tKa + eta.Ka)
        Cl <- exp(tCl)
        V <- exp(tV)
        Kemp <- exp(tKemp)
        mtime1 <- exp(tmtime1 + eta.mtime1)
        period <- tperiod
        SINE <- sin(2 * 3.14 * (time - mtime1) / period)
        EHC <- ifelse(SINE > 0, SINE, 0)
        d/dt(depot) = -Ka * depot + Kemp * GB * EHC
        d/dt(center) = Ka * depot - Cl / V * center
        d/dt(GB) = -Kemp * GB * EHC
        Cp <- center / V
        Cp ~ add(add.err) + prop(prop.err)
      })
    }
    f <- suppressMessages(ehc())
    expect_error(f$foceiModel, NA)
  })

  test_that("a failed inner evaluation is neither a gradient nor a cached value", {
    skip_on_cran()
    # sqrt(1 - eta.v) makes every prediction NaN for eta.v > 1, so likInner0()
    # fails in its observation loop, after it has reset llik and lp
    failMod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        cp <- linCmt()
        cpo <- cp * exp(eta.v) * sqrt(1 - eta.v)
        cpo ~ add(add.sd)
      })
    }
    .ui <- rxode2::assertRxUi(failMod)
    .n <- length(unique(nlmixr2data::theo_sd$ID))
    .vaeInnerSetup(.ui, nlmixr2data::theo_sd, matrix(0, .n, 1), vaeControl())
    on.exit(.vaeInnerFree(), add = TRUE)
    .good <- likInner(0.2, 1L)
    .gGood <- foceiInnerLp(0.2, 1L)
    expect_true(is.finite(.good))
    expect_true(is.finite(.gGood))
    expect_true(is.na(likInner(1.5, 1L)))
    # the failed evaluation has no gradient; its partial lp is not one
    expect_true(is.na(foceiInnerLp(1.5, 1L)))
    # and a later call at the last eta that succeeded solves again instead of
    # returning what the failed call left behind
    expect_identical(likInner(0.2, 1L), .good)
    expect_identical(foceiInnerLp(0.2, 1L), .gGood)
  })

  test_that("a log-density row the eta does not reach adds nothing to the eta gradient", {
    skip_on_cran()
    # eta.l reaches lam only through X, and X is 0 on every row of subject 3, so
    # that subject's log-density does not depend on eta.l at all.  foceiLik keeps
    # the eta-epsilon interaction branch for a non-normal endpoint.
    pois <- function() {
      ini({
        tl <- 1
        eta.l ~ 0.1
      })
      model({
        lam <- exp(tl + eta.l * X)
        y ~ pois(lam)
      })
    }
    .testSeed(42)
    .d <- do.call(
      rbind,
      lapply(1:3, function(id) {
        .t <- c(1, 2, 3, 5, 6, 8)
        .x <- if (id == 3L) rep(0, 6L) else as.numeric(.t > 4)
        data.frame(ID = id, TIME = .t, DV = stats::rpois(6L, exp(1 + 0.3 * .x)), AMT = 0, EVID = 0, X = .x)
      })
    )
    .h <- suppressWarnings(foceiLikLoad(pois, .d, "focei"))
    on.exit(foceiLikUnload(), add = TRUE)
    expect_equal(foceiLikSetThetaC_(.h$initPar), 0L)
    .eta <- matrix(c(0.1, -0.2, 0.3), .h$nid, .h$neta)
    .g <- foceiLikCondGrad_(.eta, 1L)$grad
    # each of its six rows added sqrt(DBL_EPSILON) when a zero derivative was
    # floored there
    expect_identical(.g[3, 1], 0)
    expect_identical(foceiInnerLp(0, 3L), 0)
    # the subjects the eta does reach still get the derivative of their value
    .fd <- vapply(
      1:2,
      function(i) {
        .up <- .eta
        .up[i, 1] <- .up[i, 1] + 1e-5
        .dn <- .eta
        .dn[i, 1] <- .dn[i, 1] - 1e-5
        (foceiLikCondGrad_(.up, 1L)$value[i] - foceiLikCondGrad_(.dn, 1L)$value[i]) / 2e-5
      },
      numeric(1)
    )
    expect_equal(.g[1:2, 1], .fd, tolerance = 1e-7)
  })

  test_that("a finite-differenced ETA is differenced against the prediction model", {
    skip_on_cran()
    # eventSens = "fd" finite-differences an ETA in a dosing parameter by
    # solving the prediction model at eta +/- h.  The base point was the
    # sensitivity model's: another ODE system, off by the solver error, which a
    # forward difference divides by h.  At a loose tolerance that put the
    # forward eta.f gradient of subject 1 at 11.6 against 33.8 central.
    fMod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- 0.7
        eta.f ~ 0.1
        eta.cl ~ 0.1
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        f(depot) <- exp(eta.f)
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::assertRxUi(fMod))
    .lp <- function(eventType) {
      .ctl <- .foceiInnerControl(
        list(
          rxControl = rxode2::rxControl(atol = 1e-3, rtol = 1e-3),
          sumProd = FALSE,
          optExpression = TRUE,
          literalFix = FALSE,
          addProp = "combined2",
          maxOdeRecalc = 5L,
          odeRecalcFactor = 10^0.5
        ),
        eventSens = "fd",
        eventType = eventType
      )
      .env <- .foceiInnerEnv(.ui, nlmixr2data::theo_sd, .ctl, "focei", matrix(0, 12, 2))
      foceiLikLoad_(.env)
      on.exit(foceiLikUnload_(), add = TRUE)
      vapply(1:3, function(i) foceiInnerLp(c(0.2, -0.1), i)[1], numeric(1))
    }
    expect_equal(.lp("forward"), .lp("central"), tolerance = 0.02)
  })

  test_that("an ETA reaching the prediction through lag() gets its data gradient", {
    skip_on_cran()
    # c0 is a bare symbol to symengine (lag() needs it as a real lhs); without the
    # derivative chained through lag(c0) the inner gradient holds only the prior
    # term and the inner problem pulls eta.cl to 0 whatever the data
    lagMod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- 0.7
        eta.cl ~ 0.1
        eta.f ~ 0.1
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        c0 <- central / v
        cp <- (0.5 * c0 + 0.5 * lag(c0)) * exp(eta.f)
        cp ~ add(add.sd)
      })
    }
    .ui <- rxode2::assertRxUi(lagMod)
    .n <- length(unique(nlmixr2data::theo_sd$ID))
    suppressWarnings(.vaeInnerSetup(.ui, nlmixr2data::theo_sd, matrix(0, .n, 2), vaeControl()))
    on.exit(.vaeInnerFree(), add = TRUE)
    .eta <- c(0.15, -0.1)
    .g <- foceiInnerLp(.eta, 1L)
    .fd <- vapply(
      1:2,
      function(k) {
        .h <- replace(numeric(2), k, 1e-4)
        (likInner(.eta + .h, 1L) - likInner(.eta - .h, 1L)) / 2e-4
      },
      numeric(1)
    )
    expect_equal(.g, .fd, tolerance = 1e-3)
    # fits move eta.cl off 0, and fast = TRUE declines the analytic gradient
    # (whose augmented model could not be solved: "required for solving: c0")
    # for finite differences
    for (.fast in c(FALSE, TRUE)) {
      .fit <- .nlmixr(
        lagMod,
        nlmixr2data::theo_sd,
        "focei",
        foceiControl(print = 0L, fast = .fast, maxOuterIterations = 2L, covMethod = "", calcTables = FALSE)
      )
      expect_true(is.finite(.fit$objf))
      expect_gt(stats::sd(.fit$eta$eta.cl), 0.01)
    }
  })

  test_that("a model whose only ETA reaches the prediction through lag() is fitted", {
    skip_on_cran()
    # every d(pred)/d(eta) goes through lag(c0) here, which the inner Hessian build
    # takes for "no prediction depends on a random effect" unless the derivative is
    # chained through it; the analytic covariance declines the model
    lagOnly <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- 0.7
        eta.cl ~ 0.1
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        c0 <- central / v
        cp <- 0.5 * c0 + 0.5 * lag(c0)
        cp ~ add(add.sd)
      })
    }
    .acc <- new.env(parent = emptyenv())
    .acc$msg <- character(0)
    .fit <- withCallingHandlers(
      suppressWarnings(nlmixr2(
        lagOnly,
        nlmixr2data::theo_sd,
        "focei",
        foceiControl(print = 0L, maxOuterIterations = 2L, covMethod = "analytic", calcTables = FALSE)
      )),
      message = function(m) {
        .acc$msg <- c(.acc$msg, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
    expect_true(.foceiUsesLagVar(.fit$finalUi))
    expect_true(is.finite(.fit$objf))
    expect_gt(stats::sd(.fit$eta$eta.cl), 0.1)
    expect_true(any(grepl(
      "lag() of a calculated variable is out of analytic-covariance scope",
      .acc$msg,
      fixed = TRUE
    )))
    expect_false(.covBaseName(.fit$covMethod) == "analytic")
  })

  test_that("lag() of an ODE state is refused, so only calculated variables need the lag() chain", {
    # .foceiUsesLagVar() looks at calculated variables only; if rxode2 ever accepts
    # lag() of a state, that state's sensitivities need the same chain
    .out <- utils::capture.output(
      expect_error(rxode2::rxode2("d/dt(central) <- -central\ncp <- lag(central)"), "syntax errors")
    )
    expect_true(any(grepl("state 'central': 'lag'.* not legal", .out)))
  })
})
