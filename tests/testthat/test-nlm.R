nmTest({
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

    expect_false(grepl("linear\\(wt\\)",
                       suppressMessages(rxode2::rxNorm(m$nlmRxModel$predOnly))))

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

    expect_true(grepl("linear\\(wt\\)",
                      suppressMessages(rxode2::rxNorm(m$nlmRxModel$predOnly))))
  })

  test_that("nlm makes sense", {

    dsn <- data.frame(i=1:1000)
    dsn$time <- exp(rnorm(1000))
    dsn$DV <- rbinom(1000,1,exp(-1+dsn$time)/(1+exp(-1+dsn$time)))

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

    fit2 <- .nlmixr(mod, dsn, est="nlm", nlmControl(print=0))

    expect_s3_class(fit2, "nlmixr2.nlm")

    fit2 <- .nlmixr(mod, dsn, est="bobyqa", bobyqaControl(print=0))

    expect_s3_class(fit2, "nlmixr2.bobyqa")

    fit2 <- .nlmixr(mod, dsn, est="uobyqa", uobyqaControl(print=0))

    expect_s3_class(fit2, "nlmixr2.uobyqa")

    fit2 <- .nlmixr(mod, dsn, est="newuoa", newuoaControl(print=0))

    expect_s3_class(fit2, "nlmixr2.newuoa")

    fit2 <- .nlmixr(mod, dsn, est="n1qn1", n1qn1Control(print=0))

    expect_s3_class(fit2, "nlmixr2.n1qn1")

    fit2 <- .nlmixr(mod, dsn, est="lbfgsb3c", lbfgsb3cControl(print=0))

    expect_s3_class(fit2, "nlmixr2.lbfgsb3c")

    fit3 <- suppressMessages({
      fit2 |>
        ini(g=unfix) |>
        .nlmixr(dsn, "nlm", nlmControl(solveType="grad", print=0))
    })

    expect_s3_class(fit3, "nlmixr2.nlm")

    fit4 <- suppressMessages({
      fit2 |>
        ini(g=unfix) |>
        .nlmixr(dsn, "nlm", nlmControl(solveType="fun", print=0))
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

    fit2 <- .nlmixr(one.cmt, nlmixr2data::theo_sd, est="nlm", list(print=0))

    fit1 <- .nlmixr(one.cmt, nlmixr2data::theo_sd, est="nlm",
                   nlmControl(scaleTo=0.0, scaleType="multAdd", print=0))

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

    s <- rxUiGet.nlmThetaS(list(rxode2::rxode2(mod)))
    expect_true(exists("..jacobian", envir = s, inherits = FALSE))
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
    fit <- .nlmixr(one.cmt, .dat, est="nlm", list(print=0))

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
    .sim <- rxode2::rxSolve(odeMM, .ev,
                            params = c(tka = 0.5, tvmax = log(60),
                                       tkm = log(40), tv = 3.45))
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
  })
})
