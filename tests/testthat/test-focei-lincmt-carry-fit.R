# Phase 3b.5/3b.6: fit-level validation of the linCmt() sensitivity carry
# against an independently integrated ODE reference, plus the gradient
# check across compartment configurations.  The data-aware fallbacks, the
# runtime fast path and CWRES are in test-focei-lincmt-carry-fit-fallback.R;
# shared fixtures (models, data, control) in helper-lincmt-carry.R.
#
# Slow batch (real FOCEi fits) -- see .slowBatches in tests/testthat.R.
# Everything here needs an rxode2 with the carry sentinels and skips
# cleanly on a released rxode2 without them.
#
# The ODE reference of the first test integrates under
# covsInterpolation = "nocb", as its fits are run.  linCmt() advances each
# interval with that interval's covariates under every interpolation, as an
# ODE solver does, and the carry reads them the same way; the A2 test below
# compares the two models' fits under locf, midpoint and nocb.
#
# Independence of the reference: rxSolve() on a UI auto-routes a linear
# d/dt() system back through linCmt() (useLinCmt=TRUE default), but the fit's
# inner model carries rx__sens_ states so that detection never fires (verified:
# the reference inner model has no linCmt() call, and options(rxode2.useLinCmt
# = FALSE) leaves its objective bit-identical). The DV simulation below uses a
# compiled rxode2 object (rxSolve.default, no detection) and passes
# useLinCmt=FALSE explicitly anyway.

test_that("carry fit matches the nocb linToOde ODE reference; naive does not", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  dat <- .carryFitDat()
  uiO <- rxode2::linToOde(rxode2::rxode2(.carryModCov))
  # simulate observations from the ODE truth under nocb
  m <- rxode2::rxode2("
cl = exp(tcl)*(wt/70)^0.75*exp(eta_cl)
v = exp(tv)
d/dt(central) = -(cl/v)*central
cp = central/v")
  set.seed(17)
  etaTrue <- rnorm(6, 0, 0.3)
  dv <- unlist(lapply(1:6, function(i) {
    rxode2::rxSolve(
      m,
      params = c(
        tcl = log(2),
        tv = log(20),
        eta_cl = etaTrue[i]
      ),
      events = dat[dat$id == i, ],
      returnType = "data.frame",
      covsInterpolation = "nocb",
      useLinCmt = FALSE
    )$cp
  }))
  obs <- dat$evid == 0
  dat$dv <- 0
  set.seed(99)
  dat$dv[obs] <- dv + rnorm(sum(obs), 0, 0.3)
  fit <- function(ui, carry, maxOut = 0L) {
    suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(ui, dat, est = "focei", control = .carryFitCtl(carry, maxOut))
    ))
  }
  # posthoc (fixed thetas): the pure inner-problem comparison
  fO <- fit(uiO, "none")
  fC <- fit(.carryModCov, "auto")
  fN <- fit(.carryModCov, "none")
  gapC <- abs(fC$objective - fO$objective)
  gapN <- abs(fN$objective - fO$objective)
  # measured: carry 2.5e-4 vs naive 0.145 (~580x); assert loosely
  expect_lt(gapC, 0.01)
  expect_gt(gapN, 0.05)
  expect_lt(max(abs(fC$eta$eta.cl - fO$eta$eta.cl)), 1e-3)
  # full fit: thetas/OBJF converge to the ODE reference only with the carry
  FO <- fit(uiO, "none", 200L)
  FC <- fit(.carryModCov, "auto", 200L)
  FN <- fit(.carryModCov, "none", 200L)
  expect_lt(abs(FC$objective - FO$objective), 0.01)
  expect_gt(abs(FN$objective - FO$objective), 0.05)
  expect_lt(max(abs(FC$theta - FO$theta)), 5e-3)
})

test_that("carry gradient matches FD across compartment configs; naive fails all", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  mkEv <- function(dur = 0) {
    ev <- .carryEv()
    if (dur > 0) {
      ev$rate <- ifelse(ev$evid == 1, ev$amt / dur, 0)
    }
    ev
  }
  runCfg <- function(mod, pars, ev) {
    ui <- suppressMessages(nlmixr2est::nlmixr2(mod))
    mkS <- function(carry) {
      u <- rxode2::.copyUi(ui)
      ctl <- nlmixr2est::foceiControl(linCmtSensCarry = carry)
      assign("control", ctl, envir = u)
      suppressMessages(u$foceiEnv)
    }
    sA <- mkS("auto")
    sN <- mkS("none")
    expect_true(grepl("rx_lcCarryAdv_", sA$..inner))
    mA <- suppressWarnings(rxode2::rxode2(sA$..inner))
    mN <- suppressWarnings(rxode2::rxode2(sN$..inner))
    slv <- function(mm, eta) {
      p <- pars
      p["ETA[1]"] <- eta
      rxode2::rxSolve(mm, params = p, events = ev, returnType = "data.frame")
    }
    h <- 1e-5
    fd <- (slv(mA, 0.3 + h)$rx_pred_ - slv(mA, 0.3 - h)$rx_pred_) / (2 * h)
    r0 <- slv(mA, 0.3)
    rn <- slv(mN, 0.3)
    sens <- "rx__sens_rx_pred__BY_ETA_1___"
    relC <- max(abs(r0[[sens]] - fd) / (abs(fd) + 1e-8))
    relN <- max(abs(rn[[sens]] - fd) / (abs(fd) + 1e-8))
    expect_lt(relC, 1e-6)
    expect_gt(relN, 1e-3)
  }
  # 2-cmt IV (m=2 row stride)
  runCfg(
    function() {
      ini({
      tcl <- log(2)
      tv <- log(20)
      tq <- log(1)
      tvp <- log(30)
      eta.cl ~ 0.1
      add.sd <- 0.5
    })
      model({
      cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
      v <- exp(tv)
      q <- exp(tq)
      vp <- exp(tvp)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
    },
    c(
      `THETA[1]` = log(2),
      `THETA[2]` = log(20),
      `THETA[3]` = log(1),
      `THETA[4]` = log(30),
      `THETA[5]` = 0.5,
      `ETA[1]` = 0.3
    ),
    mkEv()
  )
  # 1-cmt oral, eta+covariate on ka (slot 7, depot row)
  runCfg(
    function() {
      ini({
      tcl <- log(2)
      tv <- log(20)
      tka <- log(1.2)
      eta.ka ~ 0.1
      add.sd <- 0.5
    })
      model({
      cl <- exp(tcl)
      v <- exp(tv)
      ka <- exp(tka) * (wt / 70)^0.5 * exp(eta.ka)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
    },
    c(
      `THETA[1]` = log(2),
      `THETA[2]` = log(20),
      `THETA[3]` = log(1.2),
      `THETA[4]` = 0.5,
      `ETA[1]` = 0.3
    ),
    mkEv()
  )
  # 1-cmt IV infusion (rate history through the carry advance)
  runCfg(
    .carryModCov,
    c(
      `THETA[1]` = log(2),
      `THETA[2]` = log(20),
      `THETA[3]` = 0.5,
      `ETA[1]` = 0.3
    ),
    mkEv(dur = 2)
  )
  # 2-cmt oral (m=3)
  runCfg(
    function() {
      ini({
      tcl <- log(2)
      tv <- log(20)
      tq <- log(1)
      tvp <- log(30)
      tka <- log(1.2)
      eta.cl ~ 0.1
      add.sd <- 0.5
    })
      model({
      cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
      v <- exp(tv)
      q <- exp(tq)
      vp <- exp(tvp)
      ka <- exp(tka)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
    },
    c(
      `THETA[1]` = log(2),
      `THETA[2]` = log(20),
      `THETA[3]` = log(1),
      `THETA[4]` = log(30),
      `THETA[5]` = log(1.2),
      `THETA[6]` = 0.5,
      `ETA[1]` = 0.3
    ),
    mkEv()
  )
})

# the A2 model, oral one-compartment with cl stepping with a covariate at
# 6 h, as linCmt() and as ODEs
.carryA2Lin <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    bcv <- 0.5
    add.sd <- 0.7
    eta.cl ~ 0.1
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl + bcv * cv)
    v <- exp(tv)
    linCmt() ~ add(add.sd)
  })
}

.carryA2Ode <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    bcv <- 0.5
    add.sd <- 0.7
    eta.cl ~ 0.1
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl + bcv * cv)
    v <- exp(tv)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

.carryA2Fit <- function(mod, interp) {
  .d <- nlmixr2data::theo_sd
  .d$cv <- ifelse(.d$TIME < 6, 0, 1)
  suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    mod,
    .d,
    "focei",
    control = nlmixr2est::foceiControl(
      print = 0,
      covMethod = "",
      maxOuterIterations = 0L,
      rxControl = rxode2::rxControl(covsInterpolation = interp, atol = 1e-10, rtol = 1e-10)
    )
  )))
}

test_that("a carried linCmt() fit equals its ODE twin under locf, midpoint and nocb", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  # cl's eta is carried
  expect_identical(nrow(.foceiLinCmtCarryPairs(suppressMessages(nlmixr2est::nlmixr2(.carryA2Lin)))), 1L)
  for (.i in .carryTwinInterps()) {
    .lin <- .carryA2Fit(.carryA2Lin, .i)
    .ode <- .carryA2Fit(.carryA2Ode, .i)
    expect_equal(.lin$objf, .ode$objf, tolerance = 1e-6, label = .i)
    expect_equal(.lin$eta$eta.cl, .ode$eta$eta.cl, tolerance = 1e-5, label = .i)
    expect_equal(.lin$etaObf$OBJI, .ode$etaObf$OBJI, tolerance = 1e-6, label = .i)
  }
})

# the A2 fit by `est` at the initial estimates, with that estimator's own
# control
.carryA2FitEst <- function(mod, est, interp) {
  .d <- nlmixr2data::theo_sd
  .d$cv <- ifelse(.d$TIME < 6, 0, 1)
  .rx <- rxode2::rxControl(covsInterpolation = interp, atol = 1e-10, rtol = 1e-10)
  .ctl <- switch(
    est,
    foce = nlmixr2est::foceControl(print = 0, covMethod = "", maxOuterIterations = 0L, rxControl = .rx),
    fo = nlmixr2est::foControl(print = 0, covMethod = "", maxOuterIterations = 0L, rxControl = .rx),
    posthoc = nlmixr2est::posthocControl(print = 0, covMethod = "", rxControl = .rx)
  )
  suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(mod, .d, est, control = .ctl)))
}

test_that("foce, fo and posthoc of a carried linCmt() fit equal its ODE twin", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  # the covariate copies exist only with an rxode2 that advances an interval
  # with its own covariates
  skip_if_not(.rxFoceiLinCmtIntervalCovs(), "rxode2 advances linCmt() with the end record's covariates")
  # under locf and midpoint the carry's covariate copies precede rx_pred_ in
  # the inner model, so its rx_r_ is not column neta + 1
  for (.est in c("foce", "fo", "posthoc")) {
    for (.i in c("locf", "midpoint")) {
      .lab <- paste(.est, .i)
      .lin <- .carryA2FitEst(.carryA2Lin, .est, .i)
      .ode <- .carryA2FitEst(.carryA2Ode, .est, .i)
      expect_equal(.lin$objf, .ode$objf, tolerance = 1e-6, label = .lab)
      expect_equal(.lin$etaObf$OBJI, .ode$etaObf$OBJI, tolerance = 1e-6, label = .lab)
      expect_equal(.lin$eta$eta.cl, .ode$eta$eta.cl, tolerance = 1e-5, label = .lab)
    }
  }
})

# the covariate model declaring wt's own interpolation (nocb), as linCmt()
# and as ODEs
.carryDeclLin <- function() {
  ini({
    tcl <- log(2)
    tv <- log(20)
    eta.cl ~ 0.1
    add.sd <- 0.5
  })
  model({
    nocb(wt)
    cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
    v <- exp(tv)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

.carryDeclOde <- function() {
  ini({
    tcl <- log(2)
    tv <- log(20)
    eta.cl ~ 0.1
    add.sd <- 0.5
  })
  model({
    nocb(wt)
    cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
    v <- exp(tv)
    d / dt(central) <- -cl / v * central
    cp <- central / v
    cp ~ add(add.sd)
  })
}

# a focei fit at the initial thetas under the default (locf) interpolation
.carryDeclFit <- function(mod, data) {
  suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    mod,
    data,
    "focei",
    control = nlmixr2est::foceiControl(
      print = 0,
      maxOuterIterations = 0L,
      covMethod = "",
      rxControl = rxode2::rxControl(atol = 1e-10, rtol = 1e-10)
    )
  )))
}

test_that("a covariate's own interpolation declaration is the one the carry reads", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  .d <- .carryFitDat()
  .obs <- .d$evid == 0
  .d$dv <- 0
  .d$dv[.obs] <- 2 + 0.25 * (seq_len(sum(.obs)) %% 4)
  # wt declared nocb under the default locf
  .lin <- .carryDeclFit(.carryDeclLin, .d)
  .ode <- .carryDeclFit(.carryDeclOde, .d)
  expect_equal(.lin$objf, .ode$objf, tolerance = 1e-6)
  expect_equal(.lin$eta$eta.cl, .ode$eta$eta.cl, tolerance = 1e-5)
})

# .carryFitDat() with records that carry no wt of their own, or with wt
# changing at every record (as .carryMidEv()), and observed values
.carryMidFitDat <- function(kind) {
  .d <- .carryFitDat()
  if (kind == "infusion") {
    .d$rate <- ifelse(.d$evid == 1, .d$amt / 2, 0)
  } else if (kind == "addl") {
    .d$ii <- ifelse(.d$evid == 1, 12, 0)
    .d$addl <- ifelse(.d$evid == 1, 1L, 0L)
  } else if (kind == "na") {
    .d$wt[c(4, 12, 20)] <- NA
  } else if (kind == "every") {
    .d$wt <- 60 + 3 * stats::ave(.d$time, .d$id, FUN = seq_along)
  }
  .obs <- .d$evid == 0
  .d$dv <- 0
  .d$dv[.obs] <- 2 + 0.25 * (seq_len(sum(.obs)) %% 4)
  .d
}

# a focei fit at the initial thetas under midpoint
.carryMidFit <- function(mod, data) {
  suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    mod,
    data,
    "focei",
    control = nlmixr2est::foceiControl(
      print = 0,
      maxOuterIterations = 0L,
      covMethod = "",
      rxControl = rxode2::rxControl(covsInterpolation = "midpoint", atol = 1e-10, rtol = 1e-10)
    )
  )))
}

test_that("under midpoint a carried fit equals its ODE twin where records carry no covariate", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  skip_if_not(.rxFoceiLinCmtIntervalCovs(), "rxode2 advances linCmt() with the end record's covariates")
  .ode <- rxode2::linToOde(rxode2::rxode2(.carryModCov))
  for (.kind in c("infusion", "addl", "na", "every")) {
    .d <- .carryMidFitDat(.kind)
    .lin <- .carryMidFit(.carryModCov, .d)
    .ref <- .carryMidFit(.ode, .d)
    expect_equal(.lin$objf, .ref$objf, tolerance = 1e-6, label = .kind)
    expect_equal(.lin$eta$eta.cl, .ref$eta$eta.cl, tolerance = 1e-5, label = .kind)
  }
})

# .carryFitDat() with each subject's third dose (24 h) a reset: EVID 3 (the
# compartments reset, no dose) or EVID 4 (reset, then the dose)
.carryResetFitDat <- function(evid) {
  .d <- .carryFitDat()
  .obs <- .d$evid == 0
  .d$dv <- 0
  .d$dv[.obs] <- 2 + 0.25 * (seq_len(sum(.obs)) %% 4)
  .w <- which(.d$evid == 1 & abs(.d$time - (24 + (.d$id - 1) * 0.5)) < 1e-9)
  .d$evid[.w] <- evid
  .d
}

# .carryResetFitDat() with an observation at 22 h, after wt's step at 20 h:
# under locf the interval that ends at the reset reads the new wt
.carryResetFitDatLocf <- function(evid) {
  .d <- .carryResetFitDat(evid)
  .add <- .d[.d$evid == 0 & abs(.d$time - (15 + (.d$id - 1) * 0.5)) < 1e-9, ]
  .add$time <- .add$time + 7
  .add$wt <- .add$wt + 15
  .add$dv <- 2.5
  .d <- rbind(.d, .add)
  .d[order(.d$id, .d$time), ]
}

# a focei fit at the initial thetas under covsInterpolation `interp`
.carryResetFit <- function(mod, data, interp) {
  suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    mod,
    data,
    "focei",
    control = nlmixr2est::foceiControl(
      print = 0,
      maxOuterIterations = 0L,
      covMethod = "",
      rxControl = rxode2::rxControl(covsInterpolation = interp, atol = 1e-10, rtol = 1e-10)
    )
  )))
}

# a reset zeroes the compartments, so the carried sensitivities restart
# from 0 there rather than keep their values from before it
test_that("a carried fit with a reset equals its ODE twin", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  .carrySkipReset()
  .ode <- rxode2::linToOde(rxode2::rxode2(.carryModCov))
  for (.evid in c(3L, 4L)) {
    .d <- .carryResetFitDat(.evid)
    for (.i in c("locf", "midpoint", "nocb")) {
      .lin <- .carryResetFit(.carryModCov, .d, .i)
      .ref <- .carryResetFit(.ode, .d, .i)
      .lab <- paste("evid", .evid, .i)
      expect_equal(.lin$objf, .ref$objf, tolerance = 1e-6, label = .lab)
      expect_equal(.lin$eta$eta.cl, .ref$eta$eta.cl, tolerance = 1e-5, label = .lab)
    }
  }
})

test_that("a carried fit with a reset after a locf covariate step equals its ODE twin", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  .carrySkipReset()
  .ode <- rxode2::linToOde(rxode2::rxode2(.carryModCov))
  for (.evid in c(3L, 4L)) {
    .d <- .carryResetFitDatLocf(.evid)
    .lin <- .carryResetFit(.carryModCov, .d, "locf")
    .ref <- .carryResetFit(.ode, .d, "locf")
    .lab <- paste("evid", .evid)
    expect_equal(.lin$objf, .ref$objf, tolerance = 1e-6, label = .lab)
    expect_equal(.lin$eta$eta.cl, .ref$eta$eta.cl, tolerance = 1e-5, label = .lab)
  }
})

test_that("a carried f() and alag() fit with a reset equals its ODE twin", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  .carrySkipReset()
  .ode <- rxode2::linToOde(rxode2::rxode2(.carryModJump))
  .d <- .carryResetFitDat(3L)
  .lin <- .carryResetFit(.carryModJump, .d, "locf")
  .ref <- .carryResetFit(.ode, .d, "locf")
  expect_equal(.lin$objf, .ref$objf, tolerance = 1e-6)
  for (.e in c("eta.cl", "eta.f", "eta.lag")) {
    expect_equal(.lin$eta[[.e]], .ref$eta[[.e]], tolerance = 1e-5, label = .e)
  }
})
