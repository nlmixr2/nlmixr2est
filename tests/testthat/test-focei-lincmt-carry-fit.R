# Phase 3b.5/3b.6: fit-level validation of the linCmt() sensitivity carry.
#
# Slow batch (real FOCEi fits) -- see .slowBatches in tests/testthat.R.
# Everything here needs an rxode2 with the carry sentinels and skips
# cleanly on a released rxode2 without them.
#
# The ODE reference is generated semantics-matched: linCmt() evaluates each
# inter-row interval at the row-END covariate value (nocb), so the matched
# d/dt() reference must integrate under covsInterpolation="nocb" (verified
# against a real integrator, useLinCmt=FALSE at atol/rtol=1e-12: predictions
# agree to ~1e-11 under nocb, differ ~2.5% under locf).
#
# Independence of the reference: rxSolve() on a UI auto-routes a linear
# d/dt() system back through linCmt() (useLinCmt=TRUE default), but the fit's
# inner model carries rx__sens_ states so that detection never fires (verified:
# the reference inner model has no linCmt() call, and options(rxode2.useLinCmt
# = FALSE) leaves its objective bit-identical). The DV simulation below uses a
# compiled rxode2 object (rxSolve.default, no detection) and passes
# useLinCmt=FALSE explicitly anyway.

.carryFitDat <- function(nid = 6L) {
  do.call(rbind, lapply(seq_len(nid), function(i) {
    tim <- c(0, 3, 7, 15, 24, 30, 41, 50) + (i - 1) * 0.5
    d <- data.frame(id = i, time = tim,
                    amt = c(100, 0, 100, 0, 100, 0, 100, 0),
                    evid = c(1, 0, 1, 0, 1, 0, 1, 0), cmt = 1)
    w0 <- 60 + 5 * i
    d$wt <- ifelse(d$time < 20, w0, ifelse(d$time < 40, w0 + 15, w0 + 30))
    d
  }))
}

.carryFitCtl <- function(carry, maxOut = 0L) {
  nlmixr2est::foceiControl(print = 0, maxOuterIterations = maxOut,
                           covMethod = "", calcTables = FALSE,
                           sigdig = 8, etaNudge = 0, etaNudge2 = 0,
                           rxControl = rxode2::rxControl(covsInterpolation = "nocb"),
                           linCmtSensCarry = carry)
}

test_that("carry fit matches the nocb linToOde ODE reference; naive does not", {
  skip_on_cran()
  skip_if_not(nlmixr2est:::.rxFoceiLinCmtCarryCapable())
  fLin <- function() {
    ini({tcl <- log(2); tv <- log(20); eta.cl ~ 0.1; add.sd <- 0.5})
    model({
      cl <- exp(tcl) * (wt/70)^0.75 * exp(eta.cl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  dat <- .carryFitDat()
  uiO <- rxode2::linToOde(rxode2::rxode2(fLin))
  # simulate observations from the ODE truth under nocb
  m <- rxode2::rxode2("
cl = exp(tcl)*(wt/70)^0.75*exp(eta_cl)
v = exp(tv)
d/dt(central) = -(cl/v)*central
cp = central/v")
  set.seed(17)
  etaTrue <- rnorm(6, 0, 0.3)
  dv <- unlist(lapply(1:6, function(i) {
    rxode2::rxSolve(m, params = c(tcl = log(2), tv = log(20),
                                  eta_cl = etaTrue[i]),
                    events = dat[dat$id == i, ], returnType = "data.frame",
                    covsInterpolation = "nocb", useLinCmt = FALSE)$cp
  }))
  obs <- dat$evid == 0
  dat$dv <- 0
  set.seed(99)
  dat$dv[obs] <- dv + rnorm(sum(obs), 0, 0.3)
  fit <- function(ui, carry, maxOut = 0L) {
    suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(ui, dat, est = "focei",
                          control = .carryFitCtl(carry, maxOut))))
  }
  # posthoc (fixed thetas): the pure inner-problem comparison
  fO <- fit(uiO, "none")
  fC <- fit(fLin, "auto")
  fN <- fit(fLin, "none")
  gapC <- abs(fC$objective - fO$objective)
  gapN <- abs(fN$objective - fO$objective)
  # measured: carry 2.5e-4 vs naive 0.145 (~580x); assert loosely
  expect_lt(gapC, 0.01)
  expect_gt(gapN, 0.05)
  expect_lt(max(abs(fC$eta$eta.cl - fO$eta$eta.cl)), 1e-3)
  # full fit: thetas/OBJF converge to the ODE reference only with the carry
  FO <- fit(uiO, "none", 200L)
  FC <- fit(fLin, "auto", 200L)
  FN <- fit(fLin, "none", 200L)
  expect_lt(abs(FC$objective - FO$objective), 0.01)
  expect_gt(abs(FN$objective - FO$objective), 0.05)
  expect_lt(max(abs(FC$theta - FO$theta)), 5e-3)
})

test_that("carry gradient matches FD across compartment configs; naive fails all", {
  skip_on_cran()
  skip_if_not(nlmixr2est:::.rxFoceiLinCmtCarryCapable())
  mkEv <- function(dur = 0) {
    ev <- data.frame(id = 1, time = c(0, 3, 7, 15, 24, 30, 41, 50),
                     amt = c(100, 0, 100, 0, 100, 0, 100, 0),
                     evid = c(1, 0, 1, 0, 1, 0, 1, 0), cmt = 1)
    if (dur > 0) ev$rate <- ifelse(ev$evid == 1, ev$amt / dur, 0)
    ev$wt <- ifelse(ev$time < 20, 70, ifelse(ev$time < 40, 85, 100))
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
    sA <- mkS("auto"); sN <- mkS("none")
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
  runCfg(function() {
    ini({tcl <- log(2); tv <- log(20); tq <- log(1); tvp <- log(30)
         eta.cl ~ 0.1; add.sd <- 0.5})
    model({
      cl <- exp(tcl) * (wt/70)^0.75 * exp(eta.cl)
      v <- exp(tv); q <- exp(tq); vp <- exp(tvp)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }, c(`THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = log(1),
       `THETA[4]` = log(30), `THETA[5]` = 0.5, `ETA[1]` = 0.3), mkEv())
  # 1-cmt oral, eta+covariate on ka (slot 7, depot row)
  runCfg(function() {
    ini({tcl <- log(2); tv <- log(20); tka <- log(1.2)
         eta.ka ~ 0.1; add.sd <- 0.5})
    model({
      cl <- exp(tcl); v <- exp(tv)
      ka <- exp(tka) * (wt/70)^0.5 * exp(eta.ka)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }, c(`THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = log(1.2),
       `THETA[4]` = 0.5, `ETA[1]` = 0.3), mkEv())
  # 1-cmt IV infusion (rate history through the carry advance)
  runCfg(function() {
    ini({tcl <- log(2); tv <- log(20); eta.cl ~ 0.1; add.sd <- 0.5})
    model({
      cl <- exp(tcl) * (wt/70)^0.75 * exp(eta.cl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }, c(`THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = 0.5,
       `ETA[1]` = 0.3), mkEv(dur = 2))
  # 2-cmt oral (m=3)
  runCfg(function() {
    ini({tcl <- log(2); tv <- log(20); tq <- log(1); tvp <- log(30)
         tka <- log(1.2); eta.cl ~ 0.1; add.sd <- 0.5})
    model({
      cl <- exp(tcl) * (wt/70)^0.75 * exp(eta.cl)
      v <- exp(tv); q <- exp(tq); vp <- exp(tvp); ka <- exp(tka)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }, c(`THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = log(1),
       `THETA[4]` = log(30), `THETA[5]` = log(1.2), `THETA[6]` = 0.5,
       `ETA[1]` = 0.3), mkEv())
})

test_that("ss dose rows fall back to the standard gradient with a runInfo note", {
  skip_on_cran()
  skip_if_not(nlmixr2est:::.rxFoceiLinCmtCarryCapable())
  f <- function() {
    ini({tcl <- log(2); tv <- log(20); eta.cl ~ 0.1; add.sd <- 0.5})
    model({
      cl <- exp(tcl) * (wt/70)^0.75 * exp(eta.cl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  d <- data.frame(id = 1, time = c(0, 3, 7, 15, 24, 30),
                  amt = c(100, 0, 100, 0, 100, 0),
                  evid = c(1, 0, 1, 0, 1, 0), cmt = 1,
                  ss = c(1, 0, 0, 0, 0, 0), ii = c(12, 0, 0, 0, 0, 0))
  d$wt <- ifelse(d$time < 20, 70, 85)
  d$dv <- c(0, 3.2, 0, 2.5, 0, 2.2)
  fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(f, d, est = "focei",
                        control = .carryFitCtl("auto"))))
  expect_identical(fit$foceiControl$linCmtSensCarry, "none")
  expect_true(any(grepl("carry gradient off", unlist(fit$runInfo))))
})

test_that("runtime constant-theta fast path engages per subject and is equivalent", {
  skip_on_cran()
  skip_if_not(nlmixr2est:::.rxFoceiLinCmtCarryCapable())
  skip_if_not(exists("linCmtCarryFastStats", envir = asNamespace("rxode2"),
                     inherits = FALSE))
  .stats <- function(reset = FALSE) rxode2:::linCmtCarryFastStats(reset)
  .setFast <- function(x) rxode2:::linCmtCarrySetFast(x)
  on.exit(.setFast(TRUE), add = TRUE)
  f <- function() {
    ini({tcl <- log(2); tv <- log(20); eta.cl ~ 0.1; add.sd <- 0.5})
    model({
      cl <- exp(tcl) * (wt/70)^0.75 * exp(eta.cl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  ui <- suppressMessages(nlmixr2est::nlmixr2(f))
  u <- rxode2::.copyUi(ui)
  assign("control", .carryFitCtl("auto"), envir = u)
  s <- suppressMessages(u$foceiEnv)
  expect_true(grepl("rx_lcCarryAdv_", s$..inner))
  m <- suppressWarnings(rxode2::rxode2(s$..inner))
  pars <- c(`THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = 0.5,
            `ETA[1]` = 0.3)
  mkEv <- function(wt) {
    d <- data.frame(id = 1, time = c(0, 3, 7, 15, 24, 30, 41, 50),
                    amt = c(100, 0, 100, 0, 100, 0, 100, 0),
                    evid = c(1, 0, 1, 0, 1, 0, 1, 0), cmt = 1)
    d$wt <- wt
    d
  }
  evConst <- mkEv(70)
  evVary <- mkEv(ifelse(c(0, 3, 7, 15, 24, 30, 41, 50) < 20, 70, 85))
  .sens <- "rx__sens_rx_pred__BY_ETA_1___"
  runOne <- function(ev, fastOn) {
    .setFast(fastOn)
    invisible(.stats(TRUE))
    r <- rxode2::rxSolve(m, params = pars, events = ev,
                         returnType = "data.frame")
    list(r = r, s = .stats(FALSE))
  }
  # constant-covariate subject: every advance skips; results identical to slow
  kOn <- runOne(evConst, TRUE)
  kOff <- runOne(evConst, FALSE)
  expect_gt(kOn$s[["advCalls"]], 0)
  expect_equal(kOn$s[["advFast"]], kOn$s[["advCalls"]])
  expect_equal(kOff$s[["advFast"]], 0)
  expect_equal(kOn$r[[.sens]], kOff$r[[.sens]], tolerance = 1e-12)
  # varying subject: constant prefix skips, then permanently slow
  vOn <- runOne(evVary, TRUE)
  vOff <- runOne(evVary, FALSE)
  expect_gt(vOn$s[["advFast"]], 0)
  expect_lt(vOn$s[["advFast"]], vOn$s[["advCalls"]])
  expect_equal(vOn$r[[.sens]], vOff$r[[.sens]], tolerance = 1e-12)
})

test_that("CWRES consumes the carried gradient through the shared inner env", {
  skip_on_cran()
  skip_if_not(nlmixr2est:::.rxFoceiLinCmtCarryCapable())
  f <- function() {
    ini({tcl <- log(2); tv <- log(20); eta.cl ~ 0.1; add.sd <- 0.5})
    model({
      cl <- exp(tcl) * (wt/70)^0.75 * exp(eta.cl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  d <- .carryFitDat(3L)
  obs <- d$evid == 0
  d$dv <- 0
  d$dv[obs] <- 2 + 0.1 * seq_len(sum(obs))
  # calcTables on (default) so CWRES is computed through the shared inner env
  fitOne <- function(carry) {
    ctl <- nlmixr2est::foceiControl(print = 0, maxOuterIterations = 0L,
                                    covMethod = "", sigdig = 8,
                                    etaNudge = 0, etaNudge2 = 0,
                                    rxControl = rxode2::rxControl(covsInterpolation = "nocb"),
                                    linCmtSensCarry = carry)
    suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(f, d, est = "focei", control = ctl)))
  }
  fC <- fitOne("auto")
  fN <- fitOne("none")
  # mechanism observable: CWRES reuses the fit's own inner env (R/resid.R
  # .foceiEnv attr), and only the carry fit's inner model is carry-substituted
  innerTxt <- function(ft) paste(rxode2::rxNorm(ft$env$innerModel), collapse = "\n")
  expect_true(grepl("rx_lcCarryAdv_", innerTxt(fC)))
  expect_false(grepl("rx_lcCarryAdv_", innerTxt(fN)))
  expect_true(all(is.finite(fC$CWRES)))
  # the eta gradient feeds CWRES's linearization; the carried (correct) and
  # naive (conflated) gradients must produce different CWRES on varying-wt data
  expect_gt(max(abs(fC$CWRES - fN$CWRES)), 1e-4)
})

test_that("linear covariate interpolation on a varying eligible covariate is an error", {
  skip_on_cran()
  skip_if_not(nlmixr2est:::.rxFoceiLinCmtCarryCapable())
  f <- function() {
    ini({tcl <- log(2); tv <- log(20); eta.cl ~ 0.1; add.sd <- 0.5})
    model({
      cl <- exp(tcl) * (wt/70)^0.75 * exp(eta.cl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  d <- .carryFitDat(2L)
  d$dv <- ifelse(d$evid == 0, 3, 0)
  ctlLin <- nlmixr2est::foceiControl(
    print = 0, maxOuterIterations = 0L, covMethod = "", calcTables = FALSE,
    rxControl = rxode2::rxControl(covsInterpolation = "linear"))
  expect_error(suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(f, d, est = "focei", control = ctlLin))),
    "linear")
  # a covariate that is constant within every subject is fine under linear
  dc <- d
  dc$wt <- 70
  fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(f, dc, est = "focei", control = ctlLin)))
  expect_true(inherits(fit, "nlmixr2FitCore"))
  expect_false(grepl("rx_lcCarry", paste(rxode2::rxNorm(fit$env$innerModel), collapse = "\n")))
})

test_that("evid=2 rows fall back to the standard gradient with a runInfo note", {
  skip_on_cran()
  skip_if_not(nlmixr2est:::.rxFoceiLinCmtCarryCapable())
  f <- function() {
    ini({tcl <- log(2); tv <- log(20); eta.cl ~ 0.1; add.sd <- 0.5})
    model({
      cl <- exp(tcl) * (wt/70)^0.75 * exp(eta.cl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  d <- data.frame(id = 1, time = c(0, 3, 5, 7, 15, 24, 30),
                  amt = c(100, 0, 0, 100, 0, 100, 0),
                  evid = c(1, 0, 2, 1, 0, 1, 0), cmt = 1)
  d$wt <- ifelse(d$time < 20, 70, 85)
  d$dv <- c(0, 3.2, 0, 0, 2.5, 0, 2.2)
  fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(f, d, est = "focei",
                        control = .carryFitCtl("auto"))))
  expect_identical(fit$foceiControl$linCmtSensCarry, "none")
  expect_true(any(grepl("carry gradient off", unlist(fit$runInfo))))
  expect_false(grepl("rx_lcCarry", paste(rxode2::rxNorm(fit$env$innerModel), collapse = "\n")))
})
