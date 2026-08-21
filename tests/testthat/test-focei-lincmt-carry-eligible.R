test_that("linCmt carry eligibility detection", {

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

  # multiplicative shape, no data: candidate (varying NA)
  .p <- suppressMessages(.foceiLinCmtCarryPairs(.mult))
  expect_equal(nrow(.p), 1L)
  expect_equal(.p$slotName, "p1")
  expect_equal(.p$eta, "ETA_1_")
  expect_equal(.p$etaName, "eta.cl")
  expect_equal(.p$covs, "wt")
  expect_equal(.p$shape, "mult")
  expect_true(is.na(.p$varying))

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

  # additive shape
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

  # no covariate anywhere: not eligible
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

  # eta on a parameter without the covariate: not eligible
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

  # non-separable (eta inside the exponent): bias-to-false
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

  # same eta on two linCmt parameters: dropped (multi-slot unvalidated)
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

  # two etas in one slot's formula: dropped
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

  # second eta on a covariate-free parameter must not block the eligible one
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

  # IOV/occasion eta: dropped
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
  if (!inherits(.uiIov, "try-error")) {
    # eta.cl + iov.cl share cl's slot, so BOTH the two-eta rule and the IOV
    # rule apply; either way nothing may come back eligible
    .p <- suppressMessages(suppressWarnings(.foceiLinCmtCarryPairs(.uiIov)))
    expect_equal(nrow(.p), 0L)
  }
})
