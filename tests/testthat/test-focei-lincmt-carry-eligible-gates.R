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
  .datVary <- data.frame(id = rep(1:2, each = 4),
                         time = rep(c(0, 12, 24, 36), 2),
                         amt = rep(c(100, 0, 0, 0), 2),
                         evid = rep(c(1, 0, 0, 0), 2),
                         dv = rep(c(0, 1, 2, 1), 2),
                         wt = rep(c(70, 70, 90, 90), 2))
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
  expect_error(suppressMessages(
    .foceiLinCmtCarryPairs(.mult, data = .datVary, interpolation = "linear")),
    "linear")
  expect_error(suppressMessages(
    .foceiLinCmtCarryPairs(.mult, data = .datConst, interpolation = "linear")),
    NA)
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
      d/dt(eff) <- kin * cp - 0.5 * eff
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

test_that("a prediction that is not a bare linCmt() value is never carry-eligible", {
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
  expect_equal(nrow(.foceiLinCmtCarryPairs(ui)), 0L)
})
