nmTest({
  # Issue #286: estimation of a combined linCmt()/ODE model.  Before the
  # linCmt()->ODE translation the FOCEi dose landed in an eta-sensitivity
  # compartment, every prediction came back 0 and the objective function was
  # meaningless (~385000 instead of ~-337 at the true values).
  .lin <- function() {
    ini({
      tka <- log(1.57); tcl <- log(2.72); tv <- log(31.5)
      tke0 <- log(0.5); temax <- log(3)
      eta.ka ~ 0.2; eta.cl ~ 0.1
      pk.prop <- 0.1; pd.sd <- 0.3
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      ke0 <- exp(tke0); emax <- exp(temax)
      C2 <- linCmt()
      d/dt(ce) <- ke0 * (C2 - ce)
      eff <- emax * ce
      C2 ~ prop(pk.prop)
      eff ~ add(pd.sd)
    })
  }
  # the same model with the linear part written out; this one always worked
  .ode <- function() {
    ini({
      tka <- log(1.57); tcl <- log(2.72); tv <- log(31.5)
      tke0 <- log(0.5); temax <- log(3)
      eta.ka ~ 0.2; eta.cl ~ 0.1
      pk.prop <- 0.1; pd.sd <- 0.3
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      ke0 <- exp(tke0); emax <- exp(temax)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - (cl / v) * central
      C2 <- central / v
      d/dt(ce) <- ke0 * (C2 - ce)
      eff <- emax * ce
      C2 ~ prop(pk.prop)
      eff ~ add(pd.sd)
    })
  }

  .simData <- function(nid = 8) {
    .sim <- rxode2::rxode2({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      ke0 <- exp(tke0); emax <- exp(temax)
      C2 <- linCmt()
      d/dt(ce) <- ke0 * (C2 - ce)
      eff <- emax * ce
    })
    .th <- c(tka = log(1.57), tcl = log(2.72), tv = log(31.5),
             tke0 = log(0.5), temax = log(3))
    .ev <- rxode2::et(amt = 320, cmt = "depot")
    .ev <- rxode2::et(.ev, seq(0.5, 24, by = 2))
    .ev <- rxode2::et(.ev, id = seq_len(nid))
    .s <- rxode2::rxSolve(.sim, .th, .ev,
                          omega = lotri::lotri(eta.ka ~ 0.2, eta.cl ~ 0.1),
                          addDosing = FALSE, returnType = "data.frame")
    .obs <- rbind(
      data.frame(id = .s$id, time = .s$time,
                 dv = .s$C2 * (1 + 0.1 * stats::rnorm(nrow(.s))), cmt = "C2"),
      data.frame(id = .s$id, time = .s$time,
                 dv = .s$eff + 0.3 * stats::rnorm(nrow(.s)), cmt = "eff"))
    .obs$amt <- NA_real_
    .obs$evid <- 0
    .d <- rbind(data.frame(id = seq_len(nid), time = 0, dv = NA_real_,
                           cmt = "depot", amt = 320, evid = 1),
                .obs)
    .d[order(.d$id, .d$time, -.d$evid), ]
  }

  test_that("focei estimates a combined linCmt()/ODE model like the all-ODE model", {
    withr::with_seed(42, {
      rxode2::rxSetSeed(42)
      .d <- .simData()
    })
    .ctl <- foceiControl(print = 0, maxOuterIterations = 0, covMethod = "",
                         calcTables = FALSE)
    .fLin <- suppressWarnings(nlmixr2(.lin(), .d, est = "focei", control = .ctl))
    .fOde <- suppressWarnings(nlmixr2(.ode(), .d, est = "focei", control = .ctl))
    # the two encodings are the same statistical model, so the objective
    # function at the same (initial) values must agree
    expect_equal(.fLin$objDf$OBJF, .fOde$objDf$OBJF, tolerance = 1e-4)
  })

  # population-only versions for the nlm family (it has no etas); nlm adds a
  # theta-sensitivity state per theta, which shifted the linCmt() compartments
  # the same way the FOCEi eta-sensitivity states did
  .linPop <- function() {
    ini({
      tka <- log(1.57); tcl <- log(2.72); tv <- log(31.5)
      tke0 <- log(0.5); temax <- log(3)
      pk.prop <- 0.1; pd.sd <- 0.3
    })
    model({
      ka <- exp(tka); cl <- exp(tcl); v <- exp(tv)
      ke0 <- exp(tke0); emax <- exp(temax)
      C2 <- linCmt()
      d/dt(ce) <- ke0 * (C2 - ce)
      eff <- emax * ce
      C2 ~ prop(pk.prop)
      eff ~ add(pd.sd)
    })
  }
  .odePop <- function() {
    ini({
      tka <- log(1.57); tcl <- log(2.72); tv <- log(31.5)
      tke0 <- log(0.5); temax <- log(3)
      pk.prop <- 0.1; pd.sd <- 0.3
    })
    model({
      ka <- exp(tka); cl <- exp(tcl); v <- exp(tv)
      ke0 <- exp(tke0); emax <- exp(temax)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - (cl / v) * central
      C2 <- central / v
      d/dt(ce) <- ke0 * (C2 - ce)
      eff <- emax * ce
      C2 ~ prop(pk.prop)
      eff ~ add(pd.sd)
    })
  }

  test_that("nlm estimates a combined linCmt()/ODE model like the all-ODE model", {
    withr::with_seed(42, {
      rxode2::rxSetSeed(42)
      .d <- .simData()
    })
    # solveType="fun" uses no analytic derivatives and so never had the
    # compartment shift; "grad"/"hessian" add theta-sensitivity states and did
    for (.st in c("fun", "grad", "hessian")) {
      .ctl <- nlmControl(print = 0, iterlim = 1, calcTables = FALSE,
                         solveType = .st)
      .fLin <- suppressWarnings(nlmixr2(.linPop(), .d, est = "nlm", control = .ctl))
      .fOde <- suppressWarnings(nlmixr2(.odePop(), .d, est = "nlm", control = .ctl))
      expect_equal(.fLin$objDf$OBJF, .fOde$objDf$OBJF, tolerance = 1e-4,
                   info = paste0("solveType=", .st))
    }
  })

  test_that("the FD fallback and the table/pred models agree with the all-ODE model", {
    # The linCmt() compartments sit at a different number in the sensitivity
    # inner model than in the prediction models, which share one event table.
    # The translation removes the linCmt() block, so every model in the family
    # (inner, predNoLhs/FD fallback, predOnly/tables) uses one numbering.
    withr::with_seed(42, {
      rxode2::rxSetSeed(42)
      .d <- .simData()
    })
    # fallbackFD=TRUE routes a failed inner solve through the prediction model
    for (.fb in c(FALSE, TRUE)) {
      .ctl <- foceiControl(print = 0, maxOuterIterations = 0, covMethod = "",
                           calcTables = FALSE, fallbackFD = .fb)
      .fLin <- suppressWarnings(nlmixr2(.lin(), .d, est = "focei", control = .ctl))
      .fOde <- suppressWarnings(nlmixr2(.ode(), .d, est = "focei", control = .ctl))
      expect_equal(.fLin$objDf$OBJF, .fOde$objDf$OBJF, tolerance = 1e-4,
                   info = paste0("fallbackFD=", .fb))
    }
    # the prediction/table models are solved with the same event table
    .ctl <- foceiControl(print = 0, maxOuterIterations = 0, covMethod = "")
    .fLin <- suppressWarnings(nlmixr2(.lin(), .d, est = "focei", control = .ctl))
    .fOde <- suppressWarnings(nlmixr2(.ode(), .d, est = "focei", control = .ctl))
    for (.col in c("IPRED", "PRED", "IRES", "CWRES")) {
      expect_equal(.fLin[[.col]], .fOde[[.col]], tolerance = 1e-4, info = .col)
    }
  })

  test_that("focei and saem agree on a combined linCmt()/ODE model", {
    withr::with_seed(42, {
      rxode2::rxSetSeed(42)
      .d <- .simData(nid = 12)
    })
    .fLin <- suppressWarnings(nlmixr2(.lin(), .d, est = "focei",
                                      control = foceiControl(print = 0)))
    .fSaem <- suppressWarnings(nlmixr2(.lin(), .d, est = "saem",
                                       control = saemControl(print = 0,
                                                             nBurn = 50, nEm = 100)))
    # saem never had the compartment-numbering problem, so it is the reference
    expect_equal(unname(.fLin$parFixedDf$Estimate),
                 unname(.fSaem$parFixedDf$Estimate), tolerance = 0.15)
    # and the dose actually reaches the system now (predictions are not all 0)
    expect_gt(mean(abs(.fLin$IPRED)), 1)
  })
})
