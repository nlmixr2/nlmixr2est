# Phase 3b.5/3b.6: the linCmt() sensitivity carry through real fits --
# data-aware fallbacks (ss, linear interpolation; evid=2 rows keep the
# carry), the runtime constant-theta fast path, CWRES inheritance and
# fast=TRUE.  The reference
# comparison and the config sweep are in test-focei-lincmt-carry-fit.R;
# shared fixtures in helper-lincmt-carry.R.
#
# Slow batch (real FOCEi fits) -- see .slowBatches in tests/testthat.R.
# Everything here needs an rxode2 with the carry sentinels and skips
# cleanly on a released rxode2 without them.

.carryInnerTxt <- function(fit) {
  paste(rxode2::rxNorm(fit$env$innerModel), collapse = "\n")
}

test_that("ss dose rows fall back to the standard gradient with a runInfo note", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  d <- data.frame(
    id = 1,
    time = c(0, 3, 7, 15, 24, 30),
    amt = c(100, 0, 100, 0, 100, 0),
    evid = c(1, 0, 1, 0, 1, 0),
    cmt = 1,
    ss = c(1, 0, 0, 0, 0, 0),
    ii = c(12, 0, 0, 0, 0, 0)
  )
  d$wt <- ifelse(d$time < 20, 70, 85)
  d$dv <- c(0, 3.2, 0, 2.5, 0, 2.2)
  fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.carryModCov, d, est = "focei", control = .carryFitCtl("auto"))
  ))
  expect_identical(fit$foceiControl$linCmtSensCarry, "none")
  expect_true(any(grepl("carry gradient off", unlist(fit$runInfo))))
})

# .carryFitDat() with an EVID 2 record (a wt step at 19 h, no dose) in
# every subject
.carryEvid2Dat <- function() {
  .d <- .carryFitDat()
  .obs <- .d$evid == 0
  .d$dv <- 0
  .d$dv[.obs] <- 2 + 0.25 * (seq_len(sum(.obs)) %% 4)
  .add <- .d[.d$evid == 0 & abs(.d$time - (15 + (.d$id - 1) * 0.5)) < 1e-9, ]
  .add$time <- .add$time + 4
  .add$evid <- 2
  .add$wt <- .add$wt + 7
  .d <- rbind(.d, .add)
  .d[order(.d$id, .d$time), ]
}

# focei at the initial thetas under nocb, the one interpolation in which a
# linCmt() model and its ODE twin agree on every rxode2
.carryEvid2Fit <- function(mod, data, carry) {
  suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    mod,
    data,
    "focei",
    control = nlmixr2est::foceiControl(
      print = 0,
      maxOuterIterations = 0L,
      covMethod = "",
      linCmtSensCarry = carry,
      rxControl = rxode2::rxControl(covsInterpolation = "nocb", atol = 1e-10, rtol = 1e-10)
    )
  )))
}

test_that("evid=2 rows keep the carry, and the fit equals its ODE twin", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  .d <- .carryEvid2Dat()
  .fit <- .carryEvid2Fit(.carryModCov, .d, "auto")
  .ode <- .carryEvid2Fit(rxode2::linToOde(rxode2::rxode2(.carryModCov)), .d, "none")
  .naive <- .carryEvid2Fit(.carryModCov, .d, "none")
  expect_identical(.fit$foceiControl$linCmtSensCarry, "auto")
  expect_false(any(grepl("carry gradient off", unlist(.fit$runInfo))))
  expect_true(grepl("rx_lcCarry", .carryInnerTxt(.fit)))
  expect_equal(.fit$objf, .ode$objf, tolerance = 1e-6)
  # the standard gradient is off here, which is what the carry is for
  expect_gt(abs(.naive$objf - .ode$objf), 0.05)
})

test_that("linear covariate interpolation on a varying eligible covariate is an error", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  d <- .carryFitDat(2L)
  d$dv <- ifelse(d$evid == 0, 3, 0)
  ctlLin <- nlmixr2est::foceiControl(
    print = 0,
    maxOuterIterations = 0L,
    covMethod = "",
    calcTables = FALSE,
    rxControl = rxode2::rxControl(covsInterpolation = "linear")
  )
  expect_error(
    suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(.carryModCov, d, est = "focei", control = ctlLin)
    )),
    "linear"
  )
  # a covariate that is constant within every subject is fine under linear
  dc <- d
  dc$wt <- 70
  fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.carryModCov, dc, est = "focei", control = ctlLin)
  ))
  expect_true(inherits(fit, "nlmixr2FitCore"))
  expect_false(grepl("rx_lcCarry", .carryInnerTxt(fit)))
})

test_that("runtime constant-theta fast path engages per subject and is equivalent", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  skip_if_not(exists("linCmtCarryFastStats", envir = asNamespace("rxode2"), inherits = FALSE))
  .stats <- function(reset = FALSE) utils::getFromNamespace("linCmtCarryFastStats", "rxode2")(reset)
  .setFast <- function(x) utils::getFromNamespace("linCmtCarrySetFast", "rxode2")(x)
  on.exit(.setFast(TRUE), add = TRUE)
  u <- rxode2::.copyUi(suppressMessages(.carryUiCov()))
  assign("control", .carryFitCtl("auto"), envir = u)
  s <- suppressMessages(u$foceiEnv)
  expect_true(grepl("rx_lcCarryAdv_", s$..inner))
  m <- suppressWarnings(rxode2::rxode2(s$..inner))
  pars <- c(
    `THETA[1]` = log(2),
    `THETA[2]` = log(20),
    `THETA[3]` = 0.5,
    `ETA[1]` = 0.3
  )
  evConst <- .carryEv()
  evConst$wt <- 70
  evVary <- .carryEv()
  evVary$wt <- ifelse(evVary$time < 20, 70, 85)
  .sens <- "rx__sens_rx_pred__BY_ETA_1___"
  runOne <- function(ev, fastOn) {
    .setFast(fastOn)
    invisible(.stats(TRUE))
    r <- rxode2::rxSolve(m, params = pars, events = ev, returnType = "data.frame")
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
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  d <- .carryFitDat(3L)
  obs <- d$evid == 0
  d$dv <- 0
  d$dv[obs] <- 2 + 0.1 * seq_len(sum(obs))
  # calcTables on (default) so CWRES is computed through the shared inner env
  fitOne <- function(carry) {
    ctl <- nlmixr2est::foceiControl(
      print = 0,
      maxOuterIterations = 0L,
      covMethod = "",
      sigdig = 8,
      etaNudge = 0,
      etaNudge2 = 0,
      rxControl = rxode2::rxControl(covsInterpolation = "nocb"),
      linCmtSensCarry = carry
    )
    suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(.carryModCov, d, est = "focei", control = ctl)
    ))
  }
  fC <- fitOne("auto")
  fN <- fitOne("none")
  # mechanism observable: CWRES reuses the fit's own inner env (R/resid.R
  # .foceiEnv attr), and only the carry fit's inner model is carry-substituted
  expect_true(grepl("rx_lcCarryAdv_", .carryInnerTxt(fC)))
  expect_false(grepl("rx_lcCarryAdv_", .carryInnerTxt(fN)))
  expect_true(all(is.finite(fC$CWRES)))
  # the eta gradient feeds CWRES's linearization; the carried (correct) and
  # naive (conflated) gradients must produce different CWRES on varying-wt data
  expect_gt(max(abs(fC$CWRES - fN$CWRES)), 1e-4)
})

test_that("a carry-eligible fit survives foceiControl(fast=TRUE)", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  d <- .carryFitDat(2L)
  d$dv <- ifelse(d$evid == 0, 3, 0)
  ctl <- nlmixr2est::foceiControl(
    print = 0,
    maxOuterIterations = 2L,
    covMethod = "",
    calcTables = FALSE,
    fast = TRUE,
    rxControl = rxode2::rxControl(covsInterpolation = "nocb")
  )
  fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.carryModCov, d, est = "focei", control = ctl)
  ))
  expect_true(inherits(fit, "nlmixr2FitCore"))
  expect_true(is.finite(fit$objf))
  expect_true(grepl("rx_lcCarry", .carryInnerTxt(fit)))
})
