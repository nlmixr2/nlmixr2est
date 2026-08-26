# Phase 3b.1: carry-eligibility detection -- the model-level gates (pure
# linCmt() with a bare prediction, no direct time dependence) and the
# data-level confirmation (varying NA/TRUE/FALSE, "linear" interpolation).
# The structural shape/slot/eta rules are in
# test-focei-lincmt-carry-eligible.R.

test_that("data confirms or refutes a candidate; linear interpolation errors only when confirmed", {
  .mult <- function() {
    ini({
      tcl <- log(0.1)
      tv <- log(2)
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
  .datVary <- data.frame(
    id = rep(1:2, each = 4),
    time = rep(c(0, 12, 24, 36), 2),
    amt = rep(c(100, 0, 0, 0), 2),
    evid = rep(c(1, 0, 0, 0), 2),
    dv = rep(c(0, 1, 2, 1), 2),
    wt = rep(c(70, 70, 90, 90), 2)
  )
  .datConst <- .datVary
  .datConst$wt <- 70

  # with varying data: confirmed
  .p <- suppressMessages(.foceiLinCmtCarryPairs(.mult, data = .datVary))
  expect_equal(nrow(.p), 1L)
  expect_true(isTRUE(.p$varying))

  # constant-in-data: candidate not confirmed
  .p <- suppressMessages(.foceiLinCmtCarryPairs(.mult, data = .datConst))
  expect_equal(nrow(.p), 1L)
  expect_false(.p$varying)

  # linear interpolation: error only when confirmed varying
  expect_error(
    suppressMessages(
      .foceiLinCmtCarryPairs(.mult, data = .datVary, interpolation = "linear")
    ),
    "linear"
  )
  expect_error(
    suppressMessages(
      .foceiLinCmtCarryPairs(.mult, data = .datConst, interpolation = "linear")
    ),
    NA
  )
})

test_that("a mixed ODE + linCmt() model is never carry-eligible", {
  mixed <- function() {
    ini({
      tcl <- log(2)
      tv <- log(20)
      tkin <- log(0.5)
      eta.cl ~ 0.1
      eta.kin ~ 0.1
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
      v <- exp(tv)
      kin <- exp(tkin) * exp(eta.kin)
      cp <- linCmt()
      d / dt(eff) <- kin * cp - 0.5 * eff
      cp ~ add(add.sd)
    })
  }
  ui <- nlmixr2est::nlmixr2(mixed)
  expect_equal(nrow(.foceiLinCmtCarryPairs(ui)), 0L)
})

test_that("direct time dependence in the slot is never carry-eligible", {
  timeDep <- function() {
    ini({
      tcl <- log(2)
      tv <- log(20)
      eta.cl ~ 0.1
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl) * exp(-0.01 * t)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  ui <- nlmixr2est::nlmixr2(timeDep)
  expect_equal(nrow(.foceiLinCmtCarryPairs(ui)), 0L)
})

test_that("a prediction wrapping the linCmt() value is carried through the outer chain rule", {
  # #1004: the structural call is factored out as rx_lcConc_ and symengine
  # supplies d(pred)/d(rx_lcConc_); the carry needs the call to be unique
  scaled <- function() {
    ini({
      tcl <- log(2)
      tv <- log(20)
      eta.cl ~ 0.1
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
      v <- exp(tv)
      cp <- linCmt() * 2 + 1
      cp ~ add(add.sd)
    })
  }
  ui <- nlmixr2est::nlmixr2(scaled)
  expect_equal(nrow(.foceiLinCmtCarryPairs(ui)), 1L)
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  pars <- c(
    `THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = 0.5,
    `ETA[1]` = 0.3
  )
  r <- suppressWarnings(.carryJumpFd(scaled, pars, .carryEv(), "auto"))
  expect_lt(r$err, 1e-6)
  expect_true(grepl("rx_lcConc_~linCmtB(", r$txt, fixed = TRUE))
})

test_that("data-independent candidate pairs are memoized by the model digest", {
  .foceiLinCmtCarryMemoClear()
  .foceiLinCmtCarryMemoStats(reset = TRUE)
  .ui <- .carryUiCov()
  .p1 <- .foceiLinCmtCarryPairs(.ui)
  .st1 <- .foceiLinCmtCarryMemoStats()
  expect_equal(.st1[["misses"]], 1L)
  expect_equal(.st1[["hits"]], 0L)
  .p2 <- .foceiLinCmtCarryPairs(.ui)
  .st2 <- .foceiLinCmtCarryMemoStats()
  # mechanism: the second identical-model call is served from the memo
  expect_equal(.st2[["hits"]], 1L)
  expect_equal(.st2[["misses"]], 1L)
  expect_identical(.p1, .p2)
  # a different model gets its own digest, so a fresh detection
  .p3 <- .foceiLinCmtCarryPairs(.carryUiTwoPair())
  .st3 <- .foceiLinCmtCarryMemoStats()
  expect_equal(.st3[["misses"]], 2L)
  expect_false(identical(.p1, .p3))
  # a data-dependent call bypasses the memo entirely (no counter movement)
  .dat <- data.frame(
    id = rep(1:2, each = 3), time = rep(c(0, 12, 24), 2),
    amt = rep(c(100, 0, 0), 2), evid = rep(c(1, 0, 0), 2),
    dv = rep(c(0, 1, 2), 2), wt = rep(c(70, 90, 90), 2)
  )
  .p4 <- .foceiLinCmtCarryPairs(.ui, data = .dat)
  .st4 <- .foceiLinCmtCarryMemoStats()
  expect_equal(.st4, .st3)
  expect_true(isTRUE(.p4$varying[1]))
  # a control that keys the digest (covsInterpolation) invalidates the entry
  .ui2 <- .carryUiCov()
  rxode2::rxAssignControlValue(
    .ui2, "rxControl",
    rxode2::rxControl(covsInterpolation = "nocb")
  )
  .foceiLinCmtCarryPairs(.ui2)
  .st5 <- .foceiLinCmtCarryMemoStats()
  expect_equal(.st5[["misses"]], 3L)
})

test_that("candidate pairs persist through the rxode2 cache-directory sidecar", {
  .foceiLinCmtCarryMemoClear()
  .foceiLinCmtCarryMemoStats(reset = TRUE)
  .ui <- .carryUiCov()
  .p1 <- .foceiLinCmtCarryPairs(.ui) # compute + write the sidecar
  .key <- rxUiGet.foceiModelDigest(list(rxode2::assertRxUi(.ui)))
  .f <- .foceiLinCmtCarryCacheFile(.key)
  expect_true(file.exists(.f))
  # a "fresh session": empty session env, warm cache directory -- the next
  # call must be served from the sidecar (mechanism: fileHits moves,
  # misses does not)
  .foceiLinCmtCarryMemoClear(files = FALSE)
  .foceiLinCmtCarryMemoStats(reset = TRUE)
  .p2 <- .foceiLinCmtCarryPairs(.carryUiCov())
  .st <- .foceiLinCmtCarryMemoStats()
  expect_equal(.st[["fileHits"]], 1L)
  expect_equal(.st[["misses"]], 0L)
  expect_identical(.p1, .p2)
  # a sidecar from a different nlmixr2est build is stale: recompute + rewrite
  saveRDS(list(md5 = "not-this-build", pairs = .p1), .f)
  .foceiLinCmtCarryMemoClear(files = FALSE)
  .foceiLinCmtCarryMemoStats(reset = TRUE)
  .foceiLinCmtCarryPairs(.carryUiCov())
  .st2 <- .foceiLinCmtCarryMemoStats()
  expect_equal(.st2[["fileHits"]], 0L)
  expect_equal(.st2[["misses"]], 1L)
  expect_identical(readRDS(.f)$md5, nlmixr2.md5)
  # the opt-out bypasses the sidecar entirely
  withr::with_envvar(c(NLMIXR2EST_CARRY_MEMO = "off"), {
    .foceiLinCmtCarryMemoClear(files = FALSE)
    .foceiLinCmtCarryMemoStats(reset = TRUE)
    .foceiLinCmtCarryPairs(.carryUiCov())
    .st3 <- .foceiLinCmtCarryMemoStats()
    expect_equal(.st3[["fileHits"]], 0L)
    expect_equal(.st3[["misses"]], 0L)
  })
  .foceiLinCmtCarryMemoClear()
})
