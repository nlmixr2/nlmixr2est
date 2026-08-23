# linCmt() sensitivity carry: event-modifier (f()/alag()) channels --
# detection, generated text and the gradient
# against finite differences on eta through a real rxSolve() of the
# generated inner model.  Fit-level validation against the ODE reference is
# in test-focei-lincmt-carry-jump-fit.R; shared fixtures in
# helper-lincmt-carry.R.  Every test needs an rxode2 with the carry sentinels
# and the jump pin (which1=-8) and skips cleanly without them.

.carryJumpCapable <- function() {
  .rxFoceiLinCmtCarryCapable() && .rxFoceiLinCmtCarryJumpCapable() # nolint: object_usage_linter.
}

test_that("f()/alag() etas on a time-varying kernel become carry pairs", {
  skip_if_not(.carryJumpCapable())
  p <- .foceiLinCmtCarryPairs(.carryUiJump())
  expect_equal(nrow(p), 3L)
  # eta.cl: slot channel; eta.f: event-only (F depends on wt through
  # dlnF); eta.lag: event-only (lag on a time-varying kernel)
  expect_equal(p$slot, c(1L, NA, NA))
  expect_true(all(is.na(p$fD[-2])) && !is.na(p$fD[2]))
  expect_true(all(is.na(p$lagD[-3])) && !is.na(p$lagD[3]))
  expect_equal(p$fCmt[2], "depot")
  expect_equal(p$lagCmt[3], "depot")
})

test_that("a multiplicative f() with a covariate is exact without the carry", {
  skip_if_not(.carryJumpCapable())
  # F = h(wt) exp(eta): every input scales with exp(eta), so pred * dlnF
  # is exact across the covariate change; no pair is made for eta.f
  mod <- function() {
    ini({tcl <- log(2); tv <- log(20); tka <- log(1.2); tf <- -0.5
         eta.f ~ 0.1; add.sd <- 0.5})
    model({
      cl <- exp(tcl); v <- exp(tv); ka <- exp(tka)
      f(depot) <- exp(tf + eta.f) * (wt / 70)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  p <- .foceiLinCmtCarryPairs(suppressMessages(nlmixr2est::nlmixr2(mod)))
  expect_equal(nrow(p), 0L)
})

test_that("a lag on a constant kernel keeps the #920 term (no pair)", {
  skip_if_not(.carryJumpCapable())
  mod <- function() {
    ini({tcl <- log(2); tv <- log(20); tka <- log(1.2); tlag <- log(0.5)
         eta.lag ~ 0.1; add.sd <- 0.5})
    model({
      cl <- exp(tcl); v <- exp(tv); ka <- exp(tka)
      alag(depot) <- exp(tlag + eta.lag)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  p <- .foceiLinCmtCarryPairs(suppressMessages(nlmixr2est::nlmixr2(mod)))
  expect_equal(nrow(p), 0L)
})

test_that("a covariate-driven alag() is rejected (bias to the status quo)", {
  skip_if_not(.carryJumpCapable())
  mod <- function() {
    ini({tcl <- log(2); tv <- log(20); tka <- log(1.2); tlag <- log(0.5)
         eta.cl ~ 0.1; eta.lag ~ 0.1; add.sd <- 0.5})
    model({
      cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
      v <- exp(tv); ka <- exp(tka)
      alag(depot) <- exp(tlag + eta.lag) * (wt / 70)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  p <- .foceiLinCmtCarryPairs(suppressMessages(nlmixr2est::nlmixr2(mod)))
  expect_equal(p$eta, "ETA_1_")
})

test_that("the jump block is emitted with the pin, amounts/lag trackers and indicators", {
  skip_if_not(.carryJumpCapable())
  txt <- suppressMessages(.carryUiJump()$foceiEnv)$..inner
  expect_true(grepl(", -8, 0,", txt, fixed = TRUE))
  expect_true(grepl("rx_lcCarryPin_~", txt, fixed = TRUE))
  expect_true(grepl("rx_lcCarryD0_~(rx_lcCarryA0_-rx_lcCarryPA0_)", txt, fixed = TRUE))
  expect_true(grepl("rx_lcCarryDel_~ifelse(", txt, fixed = TRUE))
  expect_true(grepl("rx_lcCarryKP1_~", txt, fixed = TRUE))
  expect_true(grepl("rx_lcCarryUA0_~", txt, fixed = TRUE))
  expect_true(grepl("rx_lcCarryUL0_~", txt, fixed = TRUE))
  expect_equal(lengths(regmatches(txt, gregexpr("rx_lcCarryAdv_~", txt))), 1L)
  expect_true(grepl("rx__sens_rx_pred__BY_ETA_2___=rx_lcCarryS1r1_/", txt))
  expect_true(grepl("rx__sens_rx_pred__BY_ETA_3___=rx_lcCarryS2r1_/", txt))
})

test_that("a model with no jump channel emits no pin or trackers", {
  skip_if_not(.carryJumpCapable())
  txt <- suppressMessages(.carryUiCov()$foceiEnv)$..inner
  expect_false(grepl("rx_lcCarryPin_", txt, fixed = TRUE))
  expect_false(grepl("rx_lcCarryA0_", txt, fixed = TRUE))
})

test_that("f()+alag()+covariate gradients match FD for every eta; naive fails cl and lag", {
  skip_if_not(.carryJumpCapable())
  ev <- .carryEv()
  ev$cmt <- 1
  rC <- .carryJumpFd(.carryModJump, .carryJumpPars, ev)
  expect_lt(max(rC$err), 1e-6)
  rN <- .carryJumpFd(.carryModJump, .carryJumpPars, ev, carry = "none")
  expect_gt(rN$err[1], 1e-3)
  expect_gt(rN$err[3], 1e-3)
  # locf is piecewise constant too
  rL <- .carryJumpFd(.carryModJump, .carryJumpPars, ev, interp = "locf")
  expect_lt(max(rL$err), 1e-6)
})

test_that("an additive f() shape with a covariate needs the carry and matches FD", {
  skip_if_not(.carryJumpCapable())
  mod <- function() {
    ini({tcl <- log(2); tv <- log(20); tka <- log(1.2); tf <- 0
         eta.f ~ 0.1; add.sd <- 0.5})
    model({
      cl <- exp(tcl); v <- exp(tv); ka <- exp(tka)
      f(depot) <- expit(tf + eta.f + (wt - 70) / 70)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  ev <- .carryEv()
  ev$cmt <- 1
  pars <- c(`THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = log(1.2),
            `THETA[4]` = 0, `THETA[5]` = 0.5, `ETA[1]` = 0.3)
  p <- .foceiLinCmtCarryPairs(suppressMessages(nlmixr2est::nlmixr2(mod)))
  expect_equal(nrow(p), 1L)
  expect_true(isTRUE(p$fCov))
  rC <- .carryJumpFd(mod, pars, ev)
  expect_lt(max(rC$err), 1e-6)
  rN <- .carryJumpFd(mod, pars, ev, carry = "none")
  expect_gt(max(rN$err), 1e-3)
})
