# Phase 3b.3: codegen substitution for the linCmt() sensitivity carry --
# the generated model text.  Numeric gradient checks of that text are in
# test-focei-lincmt-carry-emit-fd.R; shared fixtures in helper-lincmt-carry.R.
#
# The substitution needs an rxode2 with the carry sentinels (linCmtB
# which1=-5/-6/-7, detected via linCmtCarryLiveTest) -- those tests skip
# cleanly on a released rxode2 without them.  The inert-hook test (a model
# with no eligible pair generates identical text with the carry enabled and
# disabled) runs everywhere.

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

test_that("two simultaneous pairs are both substituted with one advance", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  .inner <- .carryUiTwoPair()$foceiEnv$..inner
  expect_true(grepl("rx__sens_rx_pred__BY_ETA_1___=rx_lcCarry", .inner))
  expect_true(grepl("rx__sens_rx_pred__BY_ETA_2___=rx_lcCarry", .inner))
  expect_equal(lengths(regmatches(.inner, gregexpr("rx_lcCarryAdv_~", .inner))), 1L)
})

test_that("the naive build does not contain the carry block", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  .sN <- .carrySetControl(.carryUiCov(), "none")$foceiEnv
  expect_false(grepl("rx_lcCarryAdv_", .sN$..inner))
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
