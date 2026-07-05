nmTest({
  # matExp()/indLin() fits with default (compartment-1) dosing must match the
  # equivalent ODE model; source-first compartment order is the regression

  .mkData <- function(model, params, sd = 0.3, nid = 6, seed = 1234) {
    set.seed(seed)
    .ev <- rxode2::et(amt = 320, cmt = "depot", id = seq_len(nid)) |>
      rxode2::et(seq(0.5, 24, by = 1.5))
    .sim <- rxode2::rxSolve(model, .ev, params = params)
    .dat <- as.data.frame(.sim)[, c("id", "time", "cp")]
    .dat$cp <- .dat$cp + stats::rnorm(nrow(.dat), 0, sd)
    names(.dat) <- c("ID", "TIME", "DV")
    .dat$AMT <- 0
    .dat$EVID <- 0
    # dose row deliberately lacks cmt so it defaults to compartment 1
    .dose <- data.frame(ID = seq_len(nid), TIME = 0, DV = NA, AMT = 320, EVID = 1)
    .dat <- rbind(.dose, .dat)
    .dat[order(.dat$ID, .dat$TIME, -.dat$EVID), ]
  }

  test_that("matExp() linear model estimates identically to the ODE (focei/foce/nlm)", {
    odeLin <- function() {
      ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp ~ add(add.sd)
      })
    }
    matLin <- function() {
      ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        matExp()
        k_depot_central <- exp(tka + eta.ka)
        k_central_output <- exp(tcl) / exp(tv)
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }
    .dat <- .mkData(odeLin, c(tka = 0.6, tcl = 1.1, tv = 3.6))
    for (.met in c("focei", "foce")) {
      .fO <- .nlmixr(odeLin, .dat, est = .met, control = foceiControl(print = 0))
      .fM <- .nlmixr(matLin, .dat, est = .met, control = foceiControl(print = 0))
      expect_equal(.fM$objf, .fO$objf, tolerance = 1e-3)
      expect_equal(unname(fixef(.fM)), unname(fixef(.fO)), tolerance = 1e-3)
    }
  })

  test_that("matExp() linear model estimates identically to the ODE (agq/laplace)", {
    odeLin <- function() {
      ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp ~ add(add.sd)
      })
    }
    matLin <- function() {
      ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        matExp()
        k_depot_central <- exp(tka + eta.ka)
        k_central_output <- exp(tcl) / exp(tv)
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }
    .dat <- .mkData(odeLin, c(tka = 0.6, tcl = 1.1, tv = 3.6))
    for (.met in c("agq", "laplace")) {
      .ctl <- if (.met == "agq") agqControl(print = 0) else laplaceControl(print = 0)
      .fO <- .nlmixr(odeLin, .dat, est = .met, control = .ctl)
      .fM <- .nlmixr(matLin, .dat, est = .met, control = .ctl)
      expect_equal(.fM$objf, .fO$objf, tolerance = 1e-3)
      expect_equal(unname(fixef(.fM)), unname(fixef(.fO)), tolerance = 1e-3)
    }
  })

  test_that("matExp() + indLin() Michaelis-Menten estimates identically to the ODE (focei)", {
    odeMM <- function() {
      ini({ tka <- 0.45; tvmax <- log(60); tkm <- log(40); tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        ka <- exp(tka + eta.ka); vmax <- exp(tvmax); km <- exp(tkm); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - vmax * central / (km + central)
        cp <- central / v
        cp ~ add(add.sd)
      })
    }
    matMM <- function() {
      ini({ tka <- 0.45; tvmax <- log(60); tkm <- log(40); tv <- 3.45; eta.ka ~ 0.09; add.sd <- 0.7 })
      model({
        matExp()
        k_depot_central <- exp(tka + eta.ka)
        indLin(central) <- -exp(tvmax) * central / (exp(tkm) + central)
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }
    .dat <- .mkData(odeMM, c(tka = 0.6, tvmax = log(70), tkm = log(45), tv = 3.6))
    .fO <- .nlmixr(odeMM, .dat, est = "focei", control = foceiControl(print = 0))
    .fM <- .nlmixr(matMM, .dat, est = "focei", control = foceiControl(print = 0))
    expect_equal(.fM$objf, .fO$objf, tolerance = 1e-3)
    expect_equal(unname(fixef(.fM)), unname(fixef(.fO)), tolerance = 1e-3)
  })

  test_that("matExp() population model estimates identically to the ODE (nlm)", {
    odeMM <- function() {
      ini({ tka <- 0.45; tvmax <- log(60); tkm <- log(40); tv <- 3.45; add.sd <- 0.7 })
      model({
        ka <- exp(tka); vmax <- exp(tvmax); km <- exp(tkm); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - vmax * central / (km + central)
        cp <- central / v
        cp ~ add(add.sd)
      })
    }
    matMM <- function() {
      ini({ tka <- 0.45; tvmax <- log(60); tkm <- log(40); tv <- 3.45; add.sd <- 0.7 })
      model({
        matExp()
        k_depot_central <- exp(tka)
        indLin(central) <- -exp(tvmax) * central / (exp(tkm) + central)
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }
    .dat <- .mkData(odeMM, c(tka = 0.7, tvmax = log(80), tkm = log(55), tv = 3.6))
    .fO <- .nlmixr(odeMM, .dat, est = "nlm", list(print = 0))
    .fM <- .nlmixr(matMM, .dat, est = "nlm", list(print = 0))
    expect_equal(.fM$objf, .fO$objf, tolerance = 1e-3)
    expect_equal(unname(fixef(.fM)), unname(fixef(.fO)), tolerance = 1e-3)
  })

  test_that("matExp() estimates identically to the ODE (saem)", {
    odeS <- function() {
      ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; eta.ka ~ 0.3; add.sd <- 0.7 })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp ~ add(add.sd)
      })
    }
    matS <- function() {
      ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; eta.ka ~ 0.3; add.sd <- 0.7 })
      model({
        matExp()
        k_depot_central <- exp(tka + eta.ka)
        k_central_output <- exp(tcl) / exp(tv)
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }
    .ctl <- saemControl(print = 0, nBurn = 100, nEm = 100, seed = 42)
    .fO <- .nlmixr(odeS, nlmixr2data::theo_sd, est = "saem", control = .ctl)
    .fM <- .nlmixr(matS, nlmixr2data::theo_sd, est = "saem", control = .ctl)
    expect_equal(.fM$objf, .fO$objf, tolerance = 1e-3)
    expect_equal(unname(fixef(.fM)), unname(fixef(.fO)), tolerance = 1e-3)
  })
})
