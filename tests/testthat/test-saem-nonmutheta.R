test_that("saemControl(nonMuTheta='regress') recovers no-eta thetas and engages the mechanism", {
  skip_on_cran()

  # 1-cmt oral with a non-mu structural theta (v carries a Hill-like power with
  # no eta -> SAEM phi0 fixed effect).  tka, tv, pow all have no random effect.
  .mod <- function() {
    ini({
      tka <- log(1.2)
      tcl <- log(0.25)
      tv  <- log(5)
      pow <- 0.75
      eta.cl ~ 0.09
      prop.sd <- 0.15
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v  <- exp(tv) * (1 + pow)
      d/dt(depot)   <- -ka * depot
      d/dt(central) <-  ka * depot - cl / v * central
      cp <- central / v
      cp ~ prop(prop.sd)
    })
  }

  # simulate from truth
  .truth <- c(tka = log(1.2), tcl = log(0.25), tv = log(5), pow = 0.75)
  .times <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 24)
  .nsub <- 40
  rxode2::rxSetSeed(1042)
  .testSeed(1042)
  .et <- rxode2::et(amt = 100, time = 0, cmt = "depot")
  .et <- rxode2::et(.et, .times)
  .et <- rxode2::et(.et, id = seq_len(.nsub))
  .s <- as.data.frame(rxode2::rxSolve(.mod, .et, addDosing = FALSE))
  .obs <- data.frame(ID = .s$id, TIME = .s$time, DV = .s$sim, AMT = 0, EVID = 0)
  .dose <- data.frame(ID = seq_len(.nsub), TIME = 0, DV = NA_real_, AMT = 100, EVID = 1)
  .dat <- rbind(.dose, .obs)
  .dat <- .dat[order(.dat$ID, .dat$TIME, -.dat$EVID), ]

  .ctl <- function(mode) {
    saemControl(nBurn = 200, nEm = 150, print = 0, calcTables = FALSE,
                nonMuTheta = mode)
  }
  .fitEta <- suppressWarnings(nlmixr2(.mod, .dat, est = "saem", control = .ctl("eta")))
  .fitReg <- suppressWarnings(nlmixr2(.mod, .dat, est = "saem", control = .ctl("regress")))

  .eEta <- fixef(.fitEta)
  .eReg <- fixef(.fitReg)

  # both converge to finite estimates
  expect_true(all(is.finite(.eReg)))
  expect_true(all(is.finite(.eEta)))

  # mechanism engaged, not a silent fallback: the regress path must produce a
  # different objective than the classic phi0 path.
  expect_true(abs(as.numeric(.fitReg$objf) - as.numeric(.fitEta$objf)) > 1e-3)

  # the no-eta thetas are recovered well under regress (the target of the option)
  expect_lt(abs(.eReg[["tka"]] - .truth[["tka"]]), 0.15)
  expect_lt(abs(.eReg[["tv"]]  - .truth[["tv"]]),  0.20)
  expect_lt(abs(.eReg[["pow"]] - .truth[["pow"]]), 0.20)

  # and regress recovers the absorption theta at least as well as the phi0 path
  expect_lte(abs(.eReg[["tka"]] - .truth[["tka"]]),
             abs(.eEta[["tka"]] - .truth[["tka"]]) + 0.05)
})
