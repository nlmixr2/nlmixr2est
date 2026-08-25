# Phase 3b.1: carry-eligibility detection -- which (linCmt parameter slot,
# eta) pairs qualify structurally (shape, slot and eta rules).  The
# model-level and data-level gates are in
# test-focei-lincmt-carry-eligible-gates.R.

test_that("multiplicative shape: one eligible pair reported as a candidate", {
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
  .p <- suppressMessages(.foceiLinCmtCarryPairs(.mult))
  expect_equal(nrow(.p), 1L)
  expect_equal(.p$slotName, "p1")
  expect_equal(.p$eta, "ETA_1_")
  expect_equal(.p$etaName, "eta.cl")
  expect_equal(.p$covs, "wt")
  expect_equal(.p$shape, "mult")
  expect_true(is.na(.p$varying))
})

test_that("additive shape is eligible", {
  .add <- function() {
    ini({
      tcl <- 0.1
      tv <- log(2)
      eta.cl ~ 0.1
      add.sd <- 0.5
    })
    model({
      cl <- tcl * (wt / 70)^0.75 + eta.cl
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  .p <- suppressMessages(.foceiLinCmtCarryPairs(.add))
  expect_equal(nrow(.p), 1L)
  expect_equal(.p$shape, "add")
})

test_that("no covariate anywhere: not eligible", {
  .noCov <- function() {
    ini({
      tcl <- log(0.1)
      tv <- log(2)
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
  expect_equal(nrow(suppressMessages(.foceiLinCmtCarryPairs(.noCov))), 0L)
})

test_that("eta on a parameter without the covariate: not eligible", {
  .etaOffCov <- function() {
    ini({
      tcl <- log(0.1)
      tv <- log(2)
      eta.v ~ 0.1
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl) * (wt / 70)^0.75
      v <- exp(tv) * exp(eta.v)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  expect_equal(nrow(suppressMessages(.foceiLinCmtCarryPairs(.etaOffCov))), 0L)
})

test_that("non-separable (eta inside the exponent) biases to false", {
  .nonSep <- function() {
    ini({
      tcl <- log(0.1)
      tv <- log(2)
      eta.cl ~ 0.1
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl) * (wt / 70)^(0.75 + eta.cl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  expect_equal(nrow(suppressMessages(.foceiLinCmtCarryPairs(.nonSep))), 0L)
})

test_that("same eta on two linCmt parameters is dropped (multi-slot unvalidated)", {
  .twoSlot <- function() {
    ini({
      tcl <- log(0.1)
      tv <- log(2)
      eta.cl ~ 0.1
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
      v <- exp(tv) * exp(eta.cl)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  expect_equal(nrow(suppressMessages(.foceiLinCmtCarryPairs(.twoSlot))), 0L)
})

test_that("two etas in one slot's formula are dropped", {
  .twoEta <- function() {
    ini({
      tcl <- log(0.1)
      tv <- log(2)
      eta.cl ~ 0.1
      eta.cl2 ~ 0.1
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl + eta.cl2)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  # (two etas on one parameter also downgrades mu-referencing with a warning
  # at UI build; that is rxode2 noise, not the detection's)
  expect_equal(nrow(suppressWarnings(suppressMessages(.foceiLinCmtCarryPairs(.twoEta)))), 0L)
})

test_that("a second eta on a covariate-free parameter does not block the eligible one", {
  .mixed <- function() {
    ini({
      tcl <- log(0.1)
      tv <- log(2)
      eta.cl ~ 0.1
      eta.v ~ 0.1
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
      v <- exp(tv) * exp(eta.v)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  .p <- suppressMessages(.foceiLinCmtCarryPairs(.mixed))
  expect_equal(nrow(.p), 1L)
  expect_equal(.p$etaName, "eta.cl")
})

test_that("IOV/occasion eta is dropped", {
  .iov <- function() {
    ini({
      tcl <- log(0.1)
      tv <- log(2)
      eta.cl ~ 0.1
      iov.cl ~ 0.05 | occ
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl + iov.cl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  .uiIov <- suppressMessages(suppressWarnings(try(rxode2::rxode2(.iov), silent = TRUE)))
  skip_if(inherits(.uiIov, "try-error"))
  # eta.cl + iov.cl share cl's slot, so BOTH the two-eta rule and the IOV
  # rule apply; either way nothing may come back eligible
  .p <- suppressMessages(suppressWarnings(.foceiLinCmtCarryPairs(.uiIov)))
  expect_equal(nrow(.p), 0L)
})
