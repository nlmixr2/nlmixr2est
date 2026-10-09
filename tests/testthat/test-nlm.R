nmTest({
  test_that("nlm-family controls take only the codes there are (issue 1140)", {
    # eventType/optimHessType 1 = forward, 2 = central; solveType 1 = fun,
    # 2 = grad, 3 = hessian (optim: no hessian); there is no solve for any other
    for (.f in c("nlmControl", "nlminbControl")) {
      .ctl <- get(.f)
      expect_identical(.ctl()$eventType, 2L, info = .f)
      expect_identical(.ctl(eventType = "forward")$eventType, 1L, info = .f)
      expect_identical(.ctl(eventType = 1)$eventType, 1L, info = .f)
      expect_identical(.ctl()$optimHessType, 2L, info = .f)
      expect_identical(.ctl(optimHessType = "forward")$optimHessType, 1L, info = .f)
      expect_identical(.ctl(optimHessType = 2L)$optimHessType, 2L, info = .f)
      expect_identical(.ctl(solveType = 2)$solveType, 2L, info = .f)
      expect_identical(.ctl(solveType = "fun")$solveType, 1L, info = .f)
      for (.v in 3:6) {
        expect_error(.ctl(eventType = .v), "'eventType' must be one of", info = .f)
        expect_error(.ctl(optimHessType = .v), "'optimHessType' must be one of", info = .f)
      }
      expect_error(.ctl(solveType = 4), "'solveType' must be one of", info = .f)
      expect_error(.ctl(eventType = 1.5), "'eventType' must be one of", info = .f)
      expect_error(.ctl(eventType = NA_real_), "'eventType' must be one of", info = .f)
    }
    expect_identical(optimControl(solveType = 1)$solveType, 1L)
    expect_identical(optimControl()$solveType, 2L)
    expect_error(optimControl(solveType = 3), "'solveType' must be one of")
    expect_error(optimControl(eventType = 3), "'eventType' must be one of")
    expect_identical(optimControl(eventType = 1)$eventType, 1L)
    expect_error(nlsControl(eventType = 0), "'eventType' must be one of")
    expect_identical(nlsControl(eventType = "forward")$eventType, 1L)
    expect_error(nlmControl(eventType = "sideways"))
  })

  test_that("the nlm problem refuses codes it has no solve for (issue 1140)", {
    .mod <- function() {
      ini({
        E0 <- 0.5
        Em <- 0.5
      })
      model({
        v <- E0 + Em * time
        ll(bin) ~ DV * v - log(1 + exp(v))
      })
    }
    .d <- data.frame(ID = 1L, TIME = 1:10, AMT = 0, EVID = 0L, DV = rep(0:1, 5))
    # a control built by hand (as external engines do) skips the R checks
    for (.opt in c("eventType", "optimHessType", "solveType")) {
      .ctl <- nlmControl(print = 0L)
      .ctl[[.opt]] <- 5L
      expect_error(
        suppressMessages(nlmObjectiveSetup(.mod, .d, control = .ctl, gradient = .opt != "solveType")),
        .opt,
        info = .opt
      )
      .nlmFreeEnv()
    }
  })

  test_that("nlm models convert strings to numbers", {
    mod <- function() {
      ini({
        E0 <- 0.5
        Em <- 0.5
        E50 <- 2
        g <- fix(2)
      })
      model({
        v <- E0+Em*time^g/(E50^g+time^g)+wt
        p <- expit(v)
        if (p < 0.5) {
          a <- "good"
        }  else {
          a <- "bad"
        }
        ll(bin) ~ DV * v - log(1 + exp(v)) + (a=="good")*0.5
      })
    }

    m <- mod()

    expect_error(suppressMessages(rxode2::rxNorm(m$nlmRxModel$predOnly)), NA)
  })

  test_that("nlm models add interp", {
    mod <- function() {
      ini({
        E0 <- 0.5
        Em <- 0.5
        E50 <- 2
        g <- fix(2)
      })
      model({
        v <- E0+Em*time^g/(E50^g+time^g)+wt
        p <- expit(v)
        ll(bin) ~ DV * v - log(1 + exp(v))
      })
    }

    m <- suppressMessages(mod())

    expect_false(grepl("linear\\(wt\\)", suppressMessages(rxode2::rxNorm(m$nlmRxModel$predOnly))))

    mod <- function() {
      ini({
        E0 <- 0.5
        Em <- 0.5
        E50 <- 2
        g <- fix(2)
      })
      model({
        linear(wt)
        v <- E0+Em*time^g/(E50^g+time^g)+wt
        p <- expit(v)
        ll(bin) ~ DV * v - log(1 + exp(v))
      })
    }

    m <- suppressMessages(mod())

    expect_true(grepl("linear\\(wt\\)", suppressMessages(rxode2::rxNorm(m$nlmRxModel$predOnly))))
  })

  test_that("nlm makes sense", {
    dsn <- data.frame(i = 1:1000)
    dsn$time <- exp(rnorm(1000))
    dsn$DV <- rbinom(1000, 1, exp(-1 + dsn$time) / (1 + exp(-1 + dsn$time)))

    mod <- function() {
      ini({
        E0 <- 0.5
        Em <- 0.5
        E50 <- 2
        g <- fix(2)
      })
      model({
        v <- E0+Em*time^g/(E50^g+time^g)
        p <- expit(v)
        ll(bin) ~ DV * v - log(1 + exp(v))
      })
    }

    fit2 <- .nlmixr(mod, dsn, est = "nlm", nlmControl(print = 0))

    expect_s3_class(fit2, "nlmixr2.nlm")

    fit2 <- .nlmixr(mod, dsn, est = "bobyqa", bobyqaControl(print = 0))

    expect_s3_class(fit2, "nlmixr2.bobyqa")

    fit2 <- .nlmixr(mod, dsn, est = "uobyqa", uobyqaControl(print = 0))

    expect_s3_class(fit2, "nlmixr2.uobyqa")

    fit2 <- .nlmixr(mod, dsn, est = "newuoa", newuoaControl(print = 0))

    expect_s3_class(fit2, "nlmixr2.newuoa")

    fit2 <- .nlmixr(mod, dsn, est = "n1qn1", n1qn1Control(print = 0))

    expect_s3_class(fit2, "nlmixr2.n1qn1")

    fit2 <- .nlmixr(mod, dsn, est = "lbfgsb3c", lbfgsb3cControl(print = 0))

    expect_s3_class(fit2, "nlmixr2.lbfgsb3c")

    fit3 <- suppressMessages({
      fit2 |>
        ini(g=unfix) |>
        .nlmixr(dsn, "nlm", nlmControl(solveType = "grad", print = 0))
    })

    expect_s3_class(fit3, "nlmixr2.nlm")

    fit4 <- suppressMessages({
      fit2 |>
        ini(g=unfix) |>
        .nlmixr(dsn, "nlm", nlmControl(solveType = "fun", print = 0))
    })

    expect_s3_class(fit4, "nlmixr2.nlm")

    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- log(c(0, 2.7, 100))
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }

    fit2 <- .nlmixr(one.cmt, nlmixr2data::theo_sd, est = "nlm", list(print = 0))

    fit1 <- .nlmixr(
      one.cmt,
      nlmixr2data::theo_sd,
      est = "nlm",
      nlmControl(scaleTo = 0.0, scaleType = "multAdd", print = 0)
    )

    expect_s3_class(fit1, "nlmixr2.nlm")
  })

  test_that("matExp uses the nlm theta sensitivity path", {
    mod <- function() {
      ini({
        tka <- 0.45
        tcl <- log(c(0, 2.7, 100))
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        matExp()
        k_depot_central = exp(tka)
        k_central_output = exp(tcl) / exp(tv)
        cp = central / exp(tv)
        cp ~ add(add.sd)
      })
    }

    # A pure-linear matExp() model takes the native matrix-exponential
    # sensitivity path (#860, .sensMatExpNative()) instead of the ODE
    # flatten -- see the matching comment in test-focei-1.R.
    s <- rxUiGet.nlmThetaS(list(rxode2::rxode2(mod)))
    expect_true(isTRUE(s$..matExpNative))
    expect_false(exists("..jacobian", envir = s, inherits = FALSE))
    expect_true(exists("..sens", envir = s, inherits = FALSE))
  })

  test_that("matExp event sensitivities build HdTheta", {
    mod <- function() {
      ini({
        tka <- 0.45
        tf <- 0
        add.sd <- 0.7
      })
      model({
        matExp()
        k_depot_central = exp(tka)
        k_central_output = 0.2
        f(depot) <- expit(tf)
        cp = central / 10
        cp ~ add(add.sd)
      })
    }

    ui <- rxode2::rxode2(mod)

    s <- rxUiGet.nlmHdTheta(list(ui))
    expect_true(exists("..HdTheta", envir = s, inherits = FALSE))
    expect_true(any(grepl("rx__sens_central_BY_THETA_2__", get("..HdTheta", envir = s))))
  })

  test_that("matExp nlm model assembly runs", {
    mod <- function() {
      ini({
        tka <- 0.45
        tcl <- log(c(0, 2.7, 100))
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        matExp()
        k_depot_central = exp(tka)
        k_central_output = exp(tcl) / exp(tv)
        cp = central / exp(tv)
        cp ~ add(add.sd)
      })
    }

    env <- rxUiGet.nlmEnv(list(rxode2::rxode2(mod)))
    expect_true(exists("..nlmS", envir = env, inherits = FALSE))
    expect_true(any(grepl("k_depot_central", get("..nlmS", envir = env))))
  })

  test_that("nlm multi-subject parallel solving works", {
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- log(c(0, 2.7, 100))
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }

    .dat <- nlmixr2data::theo_md
    .nSubjects <- length(unique(.dat$ID))
    fit <- .nlmixr(one.cmt, .dat, est = "nlm", list(print = 0))

    expect_s3_class(fit, "nlmixr2.nlm")
    expect_true(.nSubjects > 1)
    expect_equal(length(unique(fit$ID)), .nSubjects)
  })

  test_that("matExp + indLin() Michaelis-Menten nlm fit matches the ODE fit", {
    # ODE Michaelis-Menten one-compartment oral model
    odeMM <- function() {
      ini({
        tka <- 0.45
        tvmax <- log(60)
        tkm <- log(40)
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        vmax <- exp(tvmax)
        km <- exp(tkm)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - vmax * central / (km + central)
        cp <- central / v
        cp ~ add(add.sd)
      })
    }
    # Equivalent matrix-exponential / inductive-linearization formulation: the
    # nonlinear Michaelis-Menten elimination is supplied via indLin().
    matMM <- function() {
      ini({
        tka <- 0.45
        tvmax <- log(60)
        tkm <- log(40)
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        matExp()
        k_depot_central <- exp(tka)
        indLin(central) <- -exp(tvmax) * central / (exp(tkm) + central)
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }

    .testSeed(123)
    .ev <- et(amt = 320, cmt = "depot", id = 1:6) |> et(seq(0.5, 24, by = 1.5))
    .sim <- rxode2::rxSolve(odeMM, .ev, params = c(tka = 0.5, tvmax = log(60), tkm = log(40), tv = 3.45))
    .dat <- as.data.frame(.sim)[, c("id", "time", "cp")]
    .dat$cp <- .dat$cp + stats::rnorm(nrow(.dat), 0, 0.3)
    names(.dat) <- c("ID", "TIME", "DV")
    .dat$AMT <- 0
    .dat$EVID <- 0
    .dose <- data.frame(ID = 1:6, TIME = 0, DV = NA, AMT = 320, EVID = 1)
    .dat <- rbind(.dose, .dat)
    .dat <- .dat[order(.dat$ID, .dat$TIME, -.dat$EVID), ]

    .fOde <- .nlmixr(odeMM, .dat, est = "nlm", list(print = 0))
    .fMat <- .nlmixr(matMM, .dat, est = "nlm", list(print = 0))

    expect_s3_class(.fMat, "nlmixr2.nlm")
    # The matExp + indLin() model and the ODE model are mathematically identical,
    # so the objective function and fixed effects must agree.
    expect_equal(.fMat$objf, .fOde$objf, tolerance = 1e-3)
    expect_equal(unname(fixef(.fMat)), unname(fixef(.fOde)), tolerance = 1e-3)
    # the ODE form is a work-around, and the fit says so
    expect_true(any(grepl("indLin() forcing", .fMat$runInfo, fixed = TRUE)))
    expect_false(any(grepl("indLin() forcing", .fOde$runInfo, fixed = TRUE)))
  })

  test_that("nlm-family covariance from a ll() model matches the Poisson GLM (#issue not doubled)", {
    # A constant-hazard RTTE ll() model is, up to a constant, the same
    # likelihood as a Poisson GLM with a log-time offset -- so their SEs
    # should agree.  Before the fix, nlm-family covariances (bobyqa/nlm/...)
    # were 4x too large (SE 2x too large) because the objective is built as
    # a plain -1*LL (not -2*LL like FOCEI/SAEM), while the Hessian->cov step
    # in `.nlmFinalizeList()` assumed a -2*LL scale.
    set.seed(1099)
    .n <- 60L
    .trt <- rep(0:1, each = .n)
    .h <- exp(-2 + -0.4 * .trt)
    .time <- pmin(stats::rexp(length(.trt), rate = .h), 5)
    .event <- as.integer(.time < 5)
    .dat <- data.frame(id = seq_along(.trt), time = .time, event = .event, trt = .trt, dv = .time, evid = 0L)

    .mod <- function() {
      ini({
        log_h0 <- log(0.2)
        log_hr <- -0.1
      })
      model({
        log_h <- log_h0 + log_hr * trt
        h <- exp(log_h)
        H <- h * time
        tte_ll <- event * log_h - H
        ll(tte) ~ tte_ll
      })
    }

    .fit <- suppressMessages(
      .nlmixr(.mod, .dat, est = "bobyqa", bobyqaControl(print = 0))
    )

    .glm <- glm(event ~ trt, offset = log(time), family = poisson(), data = .dat)

    .seNlmixr <- .fit$parFixedDf[c("log_h0", "log_hr"), "SE"]
    .seGlm <- summary(.glm)$coefficients[, "Std. Error"]

    expect_equal(unname(.seNlmixr), unname(.seGlm), tolerance = 0.1)
  })

  test_that("nlm's covMethod = \"nlm\" differences the analytic gradient (issue 1140)", {
    skip_on_cran()
    .pk <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        d / dt(depot) <- -ka * depot
        d / dt(centr) <- ka * depot - cl / v * centr
        cp <- centr / v
        cp ~ add(add.sd)
      })
    }
    # stats::nlm(hessian = TRUE) differences function values over a step of about
    # 10%, which put tcl's SE 9% low and tv's and add.sd's 9-15% high
    .own <- .nlmixr(.pk, nlmixr2data::theo_sd, est = "nlm", control = nlmControl(print = 0L, calcTables = FALSE))
    .r <- .nlmixr(.pk, nlmixr2data::theo_sd, est = "nlm", control = nlmControl(print = 0L, calcTables = FALSE, covMethod = "r"))
    expect_identical(.own$covMethod, "r (nlm)")
    expect_identical(.r$covMethod, "r")
    expect_equal(sqrt(diag(.own$cov)), sqrt(diag(.r$cov)), tolerance = 0.01)
  })
})
