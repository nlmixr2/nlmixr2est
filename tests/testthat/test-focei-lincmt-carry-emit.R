# Phase 3b.3: codegen substitution for the linCmt() sensitivity carry.
#
# The numeric tests need an rxode2 with the carry sentinels
# (linCmtB which1=-5/-6/-7, detected via linCmtCarryLiveTest) -- they skip
# cleanly on a released rxode2 without them.  The inert-hook test (a model
# with no eligible pair generates identical text with the carry enabled and
# disabled) runs everywhere.

.carryUiCov <- function() {
  one.compartment <- function() {
    ini({
      tcl <- log(2)
      tv <- log(20)
      eta.cl ~ 0.1
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  nlmixr2est::nlmixr2(one.compartment)
}

.carryUiTwoPair <- function() {
  two.pair <- function() {
    ini({
      tcl <- log(2)
      tv <- log(20)
      eta.cl ~ 0.1
      eta.v ~ 0.1
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
      v <- exp(tv) * (wt / 70) * exp(eta.v)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  nlmixr2est::nlmixr2(two.pair)
}

.carryUiNoCov <- function() {
  no.cov <- function() {
    ini({
      tcl <- log(2)
      tv <- log(20)
      eta.cl ~ 0.1
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl) * exp(eta.cl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  nlmixr2est::nlmixr2(no.cov)
}

.carrySetControl <- function(ui, value) {
  .ui <- rxode2::.copyUi(ui)
  .ctl <- nlmixr2est::foceiControl()
  .ctl$linCmtSensCarry <- value
  assign("control", .ctl, envir = .ui)
  .ui
}

# Non-uniform dose/observation spacing with a within-subject wt change
# (uniform spacing has hidden real pairing bugs in this effort before).
.carryEv <- function() {
  .ev <- data.frame(id = 1,
                    time = c(0, 3, 7, 15, 24, 30, 41, 50),
                    amt = c(100, 0, 100, 0, 100, 0, 100, 0),
                    evid = c(1, 0, 1, 0, 1, 0, 1, 0),
                    cmt = 1)
  .ev$wt <- ifelse(.ev$time < 20, 70, ifelse(.ev$time < 40, 85, 100))
  .ev
}

test_that("carry-eligible model generates the substituted carry block", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  .ui <- .carryUiCov()
  .s <- .ui$foceiEnv
  .inner <- .s$..inner
  expect_true(grepl("rx_lcCarryAdv_~linCmtB(", .inner, fixed = TRUE))
  expect_true(grepl(", -5, 0,", .inner, fixed = TRUE))
  expect_true(grepl("rx_lcCarryS0r0_~linCmtB(", .inner, fixed = TRUE))
  expect_true(grepl("rx__sens_rx_pred__BY_ETA_1___=rx_lcCarryS0r0_/", .inner))
  # the naive Jg-based line for the eligible pair must be gone
  expect_false(grepl("rx__sens_rx_pred__BY_ETA_1___=rx_expr_[0-9]+\\*linCmtB", .inner))
})

test_that("substituted gradient matches FD-on-eta; naive build shows the bug", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  .ev <- .carryEv()
  .pars <- c(`THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = 0.5,
             `ETA[1]` = 0.3)
  .solveInner <- function(inner, eta) {
    .p <- .pars
    .p["ETA[1]"] <- eta
    .m <- suppressWarnings(rxode2::rxode2(inner))
    rxode2::rxSolve(.m, params = .p, events = .ev, returnType = "data.frame")
  }
  .sC <- .carrySetControl(.carryUiCov(), "auto")$foceiEnv
  .sN <- .carrySetControl(.carryUiCov(), "none")$foceiEnv
  # the naive build must NOT contain the carry block
  expect_false(grepl("rx_lcCarryAdv_", .sN$..inner))
  .h <- 1e-5
  .rp <- .solveInner(.sC$..inner, 0.3 + .h)
  .rm <- .solveInner(.sC$..inner, 0.3 - .h)
  .fd <- (.rp$rx_pred_ - .rm$rx_pred_) / (2 * .h)
  .r0 <- .solveInner(.sC$..inner, 0.3)
  .rn <- .solveInner(.sN$..inner, 0.3)
  # predictions themselves agree between the two builds
  expect_equal(.r0$rx_pred_, .rn$rx_pred_, tolerance = 1e-10)
  .relC <- max(abs(.r0$rx__sens_rx_pred__BY_ETA_1___ - .fd) / (abs(.fd) + 1e-8))
  .relN <- max(abs(.rn$rx__sens_rx_pred__BY_ETA_1___ - .fd) / (abs(.fd) + 1e-8))
  expect_lt(.relC, 1e-6)
  # the naive gradient is measurably wrong once wt changes within the subject
  expect_gt(.relN, 1e-3)
})

test_that("two simultaneous pairs (cl and v slots) both match FD", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  .ev <- .carryEv()
  .pars <- c(`THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = 0.5,
             `ETA[1]` = 0.3, `ETA[2]` = -0.2)
  .s <- .carryUiTwoPair()$foceiEnv
  .inner <- .s$..inner
  # both etas substituted; the advance emitted exactly once
  expect_true(grepl("rx__sens_rx_pred__BY_ETA_1___=rx_lcCarry", .inner))
  expect_true(grepl("rx__sens_rx_pred__BY_ETA_2___=rx_lcCarry", .inner))
  expect_equal(lengths(regmatches(.inner, gregexpr("rx_lcCarryAdv_~", .inner))), 1L)
  .m <- suppressWarnings(rxode2::rxode2(.inner))
  .solve <- function(e1, e2) {
    .p <- .pars
    .p["ETA[1]"] <- e1
    .p["ETA[2]"] <- e2
    rxode2::rxSolve(.m, params = .p, events = .ev, returnType = "data.frame")
  }
  .h <- 1e-5
  .r0 <- .solve(0.3, -0.2)
  .fd1 <- (.solve(0.3 + .h, -0.2)$rx_pred_ - .solve(0.3 - .h, -0.2)$rx_pred_) / (2 * .h)
  .fd2 <- (.solve(0.3, -0.2 + .h)$rx_pred_ - .solve(0.3, -0.2 - .h)$rx_pred_) / (2 * .h)
  .rel1 <- max(abs(.r0$rx__sens_rx_pred__BY_ETA_1___ - .fd1) / (abs(.fd1) + 1e-8))
  .rel2 <- max(abs(.r0$rx__sens_rx_pred__BY_ETA_2___ - .fd2) / (abs(.fd2) + 1e-8))
  expect_lt(.rel1, 1e-6)
  # eta.v needs the row-local direct term on top of the amounts carry
  expect_lt(.rel2, 1e-6)
})

test_that("proportional error chains d(rx_r_)/d(eta) through the carried sensitivity", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  prop.mod <- function() {
    ini({
      tcl <- log(2)
      tv <- log(20)
      eta.cl ~ 0.1
      prop.sd <- 0.1
    })
    model({
      cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ prop(prop.sd)
    })
  }
  .s <- nlmixr2est::nlmixr2(prop.mod)$foceiEnv
  .inner <- .s$..inner
  expect_true(grepl("rx__sens_rx_r__BY_ETA_1___=.*rx__sens_rx_pred__BY_ETA_1___",
                    .inner))
  .ev <- .carryEv()
  .pars <- c(`THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = 0.1,
             `ETA[1]` = 0.3)
  .m <- suppressWarnings(rxode2::rxode2(.inner))
  .solve <- function(eta) {
    .p <- .pars
    .p["ETA[1]"] <- eta
    rxode2::rxSolve(.m, params = .p, events = .ev, returnType = "data.frame")
  }
  .h <- 1e-5
  .r0 <- .solve(0.3)
  .fdR <- (.solve(0.3 + .h)$rx_r_ - .solve(0.3 - .h)$rx_r_) / (2 * .h)
  .relR <- max(abs(.r0$rx__sens_rx_r__BY_ETA_1___ - .fdR) / (abs(.fdR) + 1e-8))
  expect_lt(.relR, 1e-6)
})

test_that("a model with no eligible pair generates identical text with carry on and off", {
  # runs regardless of rxode2 capability: with no eligible pair the hook
  # must be inert either way
  .tA <- .carrySetControl(.carryUiNoCov(), "auto")$foceiEnv$..inner
  .tN <- .carrySetControl(.carryUiNoCov(), "none")$foceiEnv$..inner
  expect_identical(.tA, .tN)
  expect_false(grepl("rx_lcCarry", .tA))
})

test_that("the model cache digest keys covsInterpolation (linear skips the carry)", {
  # runs regardless of rxode2 capability: a locf build must never be reused
  # for a linear fit of the same model, or vice versa
  .mk <- function(interp) {
    .ui <- rxode2::.copyUi(.carryUiCov())
    .ctl <- nlmixr2est::foceiControl(
      rxControl = rxode2::rxControl(covsInterpolation = interp))
    assign("control", .ctl, envir = .ui)
    .ui
  }
  expect_false(identical(.mk("locf")$foceiModelDigest,
                         .mk("linear")$foceiModelDigest))
  expect_identical(.mk("locf")$foceiModelDigest, .mk("locf")$foceiModelDigest)
})

test_that("more than four carry-eligible pairs fails loudly at model build", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  five <- function() {
    ini({
      tcl <- log(2); tv <- log(20); tq <- log(1); tv2 <- log(30)
      tq2 <- log(0.5); tv3 <- log(40); tka <- log(1)
      eta.cl ~ 0.1; eta.v ~ 0.1; eta.q ~ 0.1; eta.v2 ~ 0.1; eta.ka ~ 0.1
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
      v <- exp(tv) * (wt / 70) * exp(eta.v)
      q <- exp(tq) * (wt / 70)^0.75 * exp(eta.q)
      v2 <- exp(tv2) * (wt / 70) * exp(eta.v2)
      q2 <- exp(tq2)
      v3 <- exp(tv3)
      ka <- exp(tka) * (wt / 70) * exp(eta.ka)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  ui <- .carrySetControl(nlmixr2est::nlmixr2(five), "auto")
  expect_equal(nrow(.foceiLinCmtCarryPairs(ui)), 5L)
  s <- ui$foceiEtaS
  expect_error(.rxFoceiLinCmtCarryPairsForBuild(
    list(ui), s, paste0("ETA_", seq_len(s$..maxEta), "_")), "more than 4")
  # (inside the HdEta build the error is re-raised by the progress abort as
  # "Aborted calculation", which escapes expect_error -- the direct call
  # above is what pins the message)
})

test_that("an eta that also drives a modeled alag() keeps the #920 path, not the carry", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  lagged <- function() {
    ini({
      tka <- log(1); tcl <- log(2); tv <- log(20); tlag <- log(0.2)
      eta.cl ~ 0.1; eta.v ~ 0.1
      add.sd <- 0.5
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
      v <- exp(tv) * (wt / 70) * exp(eta.v)
      alag(depot) <- exp(tlag + eta.cl)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  ui <- .carrySetControl(nlmixr2est::nlmixr2(lagged), "auto")
  h <- ui$foceiHdEta$..HdEta
  l1 <- h[grepl("BY_ETA_1___=", h, fixed = TRUE)]
  l2 <- h[grepl("BY_ETA_2___=", h, fixed = TRUE)]
  # eta.cl (ETA_1_) drives the lag: #920 correction (-3), no carry
  expect_true(any(grepl("-3", l1, fixed = TRUE)))
  expect_false(any(grepl("rx_lcCarry", l1)))
  # eta.v (ETA_2_) is still carried
  expect_true(any(grepl("rx_lcCarry", l2)))
})
