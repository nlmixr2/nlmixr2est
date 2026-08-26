# issue 1004: the linCmt() sensitivity carry for generalized-likelihood ll()
# endpoints.  rx_pred_ is then the log-likelihood with the linCmtB() value
# call embedded in it; the carry factors the call out as rx_lcConc_,
# supplies d(rx_lcConc_)/d(eta) exactly as for a bare prediction and lets
# symengine compose the outer chain rule (and any eta dependence the
# likelihood has with the concentration held fixed).  Shared fixtures in
# helper-lincmt-carry-ll.R; fit-level checks in
# test-focei-lincmt-carry-ll-fit.R.

test_that("an explicit normal ll() endpoint is carried; the naive build is not", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  ev <- .carryLlEv()
  r <- suppressWarnings(.carryJumpFd(.carryModLlNorm, .carryLlPars, ev, "auto"))
  rn <- suppressWarnings(.carryJumpFd(.carryModLlNorm, .carryLlPars, ev, "none"))
  expect_lt(r$err, 1e-6)
  expect_gt(rn$err, 1e-3)
  l <- strsplit(r$txt, "\n")[[1]]
  # the concentration is read back once, the gradient composes the outer
  # derivative of the likelihood with the carried concentration sensitivity
  expect_equal(sum(grepl("rx_lcConc_~linCmtB(", l, fixed = TRUE)), 1L)
  expect_true(any(grepl("^rx__sens_rx_pred__BY_ETA_1___=\\(.*rx_lcConc_.*\\)\\*\\(rx_lcCarryS0r0_/", l)))
  expect_false(grepl("rx_lcConc_", rn$txt))
})

test_that("a proportional ll() (concentration used twice) is carried", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  ev <- .carryLlEv()
  pars <- .carryLlPars
  pars["THETA[3]"] <- 0.2
  r <- suppressWarnings(.carryJumpFd(.carryModLlProp, pars, ev, "auto"))
  rn <- suppressWarnings(.carryJumpFd(.carryModLlProp, pars, ev, "none"))
  expect_lt(r$err, 1e-6)
  expect_gt(rn$err, 1e-3)
  expect_equal(sum(grepl("^rx_lcConc_~", strsplit(r$txt, "\n")[[1]])), 1L)
})

test_that("a Poisson endpoint with a direct eta term outside the concentration is carried", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  ev <- .carryLlEv()
  ev$dv[ev$evid == 0] <- c(4, 3, 2, 1)
  pars <- c(`THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = 0.5,
            `THETA[4]` = 0.3, `ETA[1]` = 0.3, `ETA[2]` = -0.2)
  r <- suppressWarnings(.carryJumpFd(.carryModLlPois, pars, ev, "auto"))
  rn <- suppressWarnings(.carryJumpFd(.carryModLlPois, pars, ev, "none"))
  # eta.cl (ETA_1_) is carried; eta.b (ETA_2_) only enters outside the
  # concentration and keeps its ordinary gradient
  expect_lt(r$err[1], 1e-6)
  expect_lt(r$err[2], 1e-6)
  expect_gt(rn$err[1], 1e-3)
  l <- strsplit(r$txt, "\n")[[1]]
  expect_true(any(grepl("^rx__sens_rx_pred__BY_ETA_1___=\\(", l)))
  expect_false(any(grepl("^rx__sens_rx_pred__BY_ETA_2___=.*rx_lcCarry", l)))
})

test_that("an eta both in the kernel and outside the concentration gets both terms", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  ev <- .carryLlEv()
  ev$dv[ev$evid == 0] <- c(4, 3, 2, 1)
  pars <- c(`THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = 0.5,
            `THETA[4]` = 0.3, `ETA[1]` = 0.3)
  r <- suppressWarnings(.carryJumpFd(.carryModLlPoisShared, pars, ev, "auto"))
  expect_lt(r$err, 1e-6)
  l <- strsplit(r$txt, "\n")[[1]]
  # the composed line carries the direct term as a separate summand
  expect_true(any(grepl("^rx__sens_rx_pred__BY_ETA_1___=\\(.*\\)\\*\\(rx_lcCarryS0r0_/.*\\)\\+\\(", l)))
})

test_that("an ll() model without a time-varying covariate generates identical text", {
  .tA <- .carrySetControl(suppressMessages(nlmixr2est::nlmixr2(.carryModLlNoCov)), "auto")$foceiEnv$..inner
  .tN <- .carrySetControl(suppressMessages(nlmixr2est::nlmixr2(.carryModLlNoCov)), "none")$foceiEnv$..inner
  expect_identical(.tA, .tN)
  expect_false(grepl("rx_lcConc_", .tA))
})

test_that("the structural call must be unique and in value form", {
  mk <- function(txt) {
    s <- new.env()
    assign("rx_pred_", symengine::S(txt), envir = s) # nolint: object_name_linter.
    s
  }
  one <- "linCmtB(rx__PTR__, t, 1.0, 1.0, 0.0, -1.0, -1.0, 1.0, a, b, 0.0, 0.0, 0.0, 0.0, 0.0)"
  two <- "linCmtB(rx__PTR__, t, 1.0, 1.0, 0.0, -1.0, -1.0, 1.0, a, c, 0.0, 0.0, 0.0, 0.0, 0.0)"
  sens <- "linCmtB(rx__PTR__, t, 1.0, 1.0, 0.0, -2.0, 0.0, 1.0, a, b, 0.0, 0.0, 0.0, 0.0, 0.0)"
  pc <- .rxFoceiCarryPredCall(mk(paste0("log(", one, ") + 2*", one)))
  expect_false(pc$bare)
  expect_equal(paste(pc$predSym), "2*rx_lcConc_ + log(rx_lcConc_)")
  expect_true(.rxFoceiCarryPredCall(mk(one))$bare)
  expect_null(.rxFoceiCarryPredCall(mk(paste0(one, " + ", two))))
  expect_null(.rxFoceiCarryPredCall(mk(sens)))
  expect_null(.rxFoceiCarryPredCall(mk("a + b")))
})

test_that("a carried ll() model keeps the finite-difference inner Hessian under fast=TRUE", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  ui <- rxode2::.copyUi(suppressMessages(nlmixr2est::nlmixr2(.carryModLlNorm)))
  assign("control", nlmixr2est::foceiControl(fast = TRUE, linCmtSensCarry = "auto"),
         envir = ui)
  s <- suppressWarnings(suppressMessages(ui$foceEnv))
  expect_false(is.null(s$..linCmtCarryPairs))
  expect_null(s$..HdEta2)
})
