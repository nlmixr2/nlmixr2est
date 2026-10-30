# The data-aware carry gate for records that modify a linCmt() compartment
# outright (R/foceiLinCmtEventGate.R): which records the data has, why the
# carry cannot run on them with the loaded rxode2lincmt, and when the
# standard gradient is approximate too.  Shared fixtures in
# helper-lincmt-carry.R.

# stand-ins for the probe mocks
.gateTrue <- function() TRUE
.gateFalse <- function() FALSE

.gateKinds <- function(evid) unlist(.rxFoceiLinCmtEventKinds(evid))

# a dose, an observation, `evid`, an observation; `ss` on the dose
.gateRd <- function(evid, ss = 0) data.frame(EVID = c(1, 0, evid, 0), SS = c(ss, 0, 0, 0))

test_that("the event kinds come from standard and encoded EVIDs", {
  expect_identical(.gateKinds(c(0, 1, 2, 7)), c(reset = FALSE, repMult = FALSE))
  expect_identical(.gateKinds(c(0, 1, 3)), c(reset = TRUE, repMult = FALSE))
  expect_identical(.gateKinds(c(0, 4)), c(reset = TRUE, repMult = FALSE))
  expect_identical(.gateKinds(c(0, 1, 5)), c(reset = FALSE, repMult = TRUE))
  expect_identical(.gateKinds(c(0, 1, 6)), c(reset = FALSE, repMult = TRUE))
  # rxode2's encoded EVIDs: a bolus into cmt 1, a replace and a multiply
  expect_identical(.gateKinds(c(0, 101)), c(reset = FALSE, repMult = FALSE))
  expect_identical(.gateKinds(c(0, 40101)), c(reset = FALSE, repMult = TRUE))
  expect_identical(.gateKinds(c(0, 50101)), c(reset = FALSE, repMult = TRUE))
  expect_identical(.gateKinds(c(NA, 0)), c(reset = FALSE, repMult = FALSE))
  expect_identical(.gateKinds(NULL), c(reset = FALSE, repMult = FALSE))
})

test_that("the gate names why the carry cannot run, case by case", {
  local_mocked_bindings(
    .rxFoceiLinCmtCarryResetOk = .gateFalse,
    .rxFoceiLinCmtCarryEventOk = .gateFalse,
    .rxFoceiLinCmtJacEventOk = .gateFalse
  )
  expect_identical(.rxFoceiLinCmtCarryEventWhy(.gateRd(0, ss = 1)), "ss rows")
  expect_null(.rxFoceiLinCmtCarryEventWhy(.gateRd(0)))
  # evid=2 rows keep the carry
  expect_null(.rxFoceiLinCmtCarryEventWhy(.gateRd(2)))
  expect_identical(.rxFoceiLinCmtCarryEventWhy(.gateRd(3)), "reset rows")
  expect_identical(.rxFoceiLinCmtCarryEventWhy(.gateRd(4)), "reset rows")
  expect_identical(.rxFoceiLinCmtCarryEventWhy(.gateRd(5)), "replace/multiply rows")
  expect_identical(.rxFoceiLinCmtCarryEventWhy(.gateRd(6)), "replace/multiply rows")
  local_mocked_bindings(.rxFoceiLinCmtCarryResetOk = .gateTrue)
  expect_null(.rxFoceiLinCmtCarryEventWhy(.gateRd(3)))
  # a replace needs both the carry's and the Jacobian's handling
  local_mocked_bindings(.rxFoceiLinCmtCarryEventOk = .gateTrue)
  expect_identical(.rxFoceiLinCmtCarryEventWhy(.gateRd(5)), "replace/multiply rows")
  local_mocked_bindings(.rxFoceiLinCmtCarryEventOk = .gateFalse, .rxFoceiLinCmtJacEventOk = .gateTrue)
  expect_identical(.rxFoceiLinCmtCarryEventWhy(.gateRd(5)), "replace/multiply rows")
  local_mocked_bindings(.rxFoceiLinCmtCarryEventOk = .gateTrue)
  expect_null(.rxFoceiLinCmtCarryEventWhy(.gateRd(5)))
  expect_null(.rxFoceiLinCmtCarryEventWhy(.gateRd(6)))
})

test_that("a linCmt() fit's eta gradients are approximate only with replace/multiply and no Jacobian fix", {
  .lin <- suppressMessages(nlmixr2est::nlmixr2(.carryModCov))
  .ode <- suppressMessages(nlmixr2est::nlmixr2(rxode2::linToOde(rxode2::rxode2(.carryModCov))))
  local_mocked_bindings(.rxFoceiLinCmtJacEventOk = .gateFalse)
  expect_true(.rxFoceiLinCmtJacApprox(.lin, .gateRd(5)))
  expect_true(.rxFoceiLinCmtJacApprox(.lin, .gateRd(6)))
  expect_false(.rxFoceiLinCmtJacApprox(.lin, .gateRd(3)))
  expect_false(.rxFoceiLinCmtJacApprox(.ode, .gateRd(5)))
  local_mocked_bindings(.rxFoceiLinCmtJacEventOk = .gateTrue)
  expect_false(.rxFoceiLinCmtJacApprox(.lin, .gateRd(5)))
})

test_that("the event probes answer for the loaded rxode2lincmt, once", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  for (.f in list(.rxFoceiLinCmtCarryResetOk, .rxFoceiLinCmtCarryEventOk, .rxFoceiLinCmtJacEventOk)) {
    .a <- .f()
    expect_true(is.logical(.a) && length(.a) == 1L && !is.na(.a))
    expect_identical(.f(), .a)
  }
  # the carry probes leave the fast-path setting as it was
  .old <- rxode2::linCmtCarrySetFast(FALSE)
  on.exit(rxode2::linCmtCarrySetFast(.old), add = TRUE)
  rm(list = ls(.linCmtEventEnv), envir = .linCmtEventEnv)
  .rxFoceiLinCmtCarryResetOk()
  .rxFoceiLinCmtCarryEventOk()
  expect_false(rxode2::linCmtCarrySetFast(FALSE))
})

# .carryFitDat() with each subject's 24 h dose made `evid` (amount `amt`)
.gateEventDat <- function(evid, amt) {
  .d <- .carryFitDat()
  .obs <- .d$evid == 0
  .d$dv <- 0
  .d$dv[.obs] <- 2 + 0.25 * (seq_len(sum(.obs)) %% 4)
  .w <- which(.d$evid == 1 & abs(.d$time - (24 + (.d$id - 1) * 0.5)) < 1e-9)
  .d$evid[.w] <- evid
  .d$amt[.w] <- amt
  .d
}

# focei at the initial thetas under nocb, where a linCmt() model and its ODE
# twin agree on every rxode2
.gateFit <- function(mod, data, carry = "auto") {
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

.gateInner <- function(fit) paste(rxode2::rxNorm(fit$env$innerModel), collapse = "\n")

test_that("a reset keeps the carry where it restarts, and equals the ODE twin", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  skip_if_not(.rxFoceiLinCmtCarryResetOk(), "the loaded rxode2lincmt's carry does not restart at a reset")
  .ode <- rxode2::linToOde(rxode2::rxode2(.carryModCov))
  for (.e in c(3L, 4L)) {
    .d <- .gateEventDat(.e, if (.e == 3L) 0 else 100)
    .fit <- .gateFit(.carryModCov, .d)
    expect_true(grepl("rx_lcCarry", .gateInner(.fit)), label = .e)
    expect_equal(.fit$objf, .gateFit(.ode, .d, "none")$objf, tolerance = 1e-6, label = .e)
  }
})

test_that("a reset turns the carry off where it does not restart, and says so", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  skip_if(.rxFoceiLinCmtCarryResetOk(), "the loaded rxode2lincmt's carry restarts at a reset")
  .fit <- .gateFit(.carryModCov, .gateEventDat(3L, 0))
  expect_identical(.fit$foceiControl$linCmtSensCarry, "none")
  expect_true(any(unlist(.fit$runInfo) == "reset rows: linCmt() carry gradient off for this fit"))
})

test_that("replace and multiply keep the carry where both fixes are in, and equal the ODE twin", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  skip_if_not(
    .rxFoceiLinCmtCarryEventOk() && .rxFoceiLinCmtJacEventOk(),
    "the loaded rxode2lincmt does not follow replace/multiply records"
  )
  .ode <- rxode2::linToOde(rxode2::rxode2(.carryModCov))
  for (.e in c(5L, 6L)) {
    .d <- .gateEventDat(.e, if (.e == 5L) 50 else 0.5)
    .fit <- .gateFit(.carryModCov, .d)
    expect_false(any(grepl("replace/multiply", unlist(.fit$runInfo))), label = .e)
    expect_true(grepl("rx_lcCarry", .gateInner(.fit)), label = .e)
    expect_equal(.fit$objf, .gateFit(.ode, .d, "none")$objf, tolerance = 1e-6, label = .e)
  }
})

test_that("replace and multiply say the gradient is approximate without the Jacobian fix", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  skip_if(.rxFoceiLinCmtJacEventOk(), "the loaded rxode2lincmt's Jacobian follows replace/multiply")
  .fit <- .gateFit(.carryModCov, .gateEventDat(5L, 50))
  .ri <- unlist(.fit$runInfo)
  expect_identical(.fit$foceiControl$linCmtSensCarry, "none")
  expect_true(any(.ri == "replace/multiply rows: linCmt() carry gradient off for this fit"))
  expect_true(any(.ri == "replace/multiply rows: linCmt() eta gradients are approximate"))
  # also with the carry off by request
  .ri <- unlist(.gateFit(.carryModCov, .gateEventDat(6L, 0.5), "none")$runInfo)
  expect_true(any(.ri == "replace/multiply rows: linCmt() eta gradients are approximate"))
})
