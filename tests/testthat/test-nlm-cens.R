nmTest({
  # Tests for censoring support in the NLM estimation engine
  # These tests verify that NLM handles M2, M3, M4 censoring
  # matching the behavior used in FOCEI/SAEM

  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- log(c(0, 2.7, 100))
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  # Base dataset (no censoring)
  .dat <- nlmixr2data::theo_sd

  test_that("nlm works without censoring (no CENS column)", {
    fit <- suppressMessages(suppressWarnings(
      .nlmixr(one.cmt, .dat, est = "nlm", list(print = 0))
    ))
    expect_s3_class(fit, "nlmixr2.nlm")
    expect_equal(as.character(fit$censInformation), "No censoring")
  })

  test_that("nlm works with CENS column all zeros", {
    .dat0 <- .dat
    .dat0$CENS <- 0L
    fit0 <- suppressMessages(suppressWarnings(
      .nlmixr(one.cmt, .dat0, est = "nlm", list(print = 0))
    ))
    expect_s3_class(fit0, "nlmixr2.nlm")
    expect_equal(as.character(fit0$censInformation), "No censoring")
  })

  # Create censored datasets
  # Mark observations below a threshold as left-censored (M3)
  .LLOQ <- 2.0
  .datM3 <- .dat
  .datM3$CENS <- ifelse(.datM3$DV < .LLOQ & .datM3$EVID == 0, 1L, 0L)
  # For censored obs, set DV to LLOQ
  .datM3$DV[.datM3$CENS == 1] <- .LLOQ

  test_that("nlm accepts and processes M3 (left) censored data", {
    fit_m3 <- suppressMessages(suppressWarnings(
      .nlmixr(one.cmt, .datM3, est = "nlm", list(print = 0))
    ))
    expect_s3_class(fit_m3, "nlmixr2.nlm")
    expect_equal(as.character(fit_m3$censInformation), "M3 censoring")
  })

  test_that("nlm M3 censoring changes the objective function vs no censoring", {
    fit_base <- suppressMessages(suppressWarnings(
      .nlmixr(one.cmt, .dat, est = "nlm", list(print = 0))
    ))
    fit_m3 <- suppressMessages(suppressWarnings(
      .nlmixr(one.cmt, .datM3, est = "nlm", list(print = 0))
    ))
    expect_false(isTRUE(all.equal(fit_base$objf, fit_m3$objf)))
  })

  # M2 censoring: CENS=0 with a finite LIMIT
  .datM2 <- .dat
  .datM2$CENS <- 0L
  .datM2$LIMIT <- 0  # interval censoring: all obs have a lower bound of 0

  test_that("nlm accepts and processes M2 (interval) censored data", {
    fit_m2 <- suppressMessages(suppressWarnings(
      .nlmixr(one.cmt, .datM2, est = "nlm", list(print = 0))
    ))
    expect_s3_class(fit_m2, "nlmixr2.nlm")
    expect_equal(as.character(fit_m2$censInformation), "M2 censoring")
  })

  test_that("nlm M2 censoring changes the objective function vs no censoring", {
    fit_base <- suppressMessages(suppressWarnings(
      .nlmixr(one.cmt, .dat, est = "nlm", list(print = 0))
    ))
    fit_m2 <- suppressMessages(suppressWarnings(
      .nlmixr(one.cmt, .datM2, est = "nlm", list(print = 0))
    ))
    expect_false(isTRUE(all.equal(fit_base$objf, fit_m2$objf)))
  })

  # M4 censoring: CENS!=0 with a finite LIMIT
  .datM4 <- .datM3
  .datM4$LIMIT <- 0  # add LIMIT for M4

  test_that("nlm accepts and processes M4 (interval-censored) data", {
    fit_m4 <- suppressMessages(suppressWarnings(
      .nlmixr(one.cmt, .datM4, est = "nlm", list(print = 0))
    ))
    expect_s3_class(fit_m4, "nlmixr2.nlm")
    expect_equal(as.character(fit_m4$censInformation), "M2 and M4 censoring")
  })

  test_that("nlm M3 and M4 give different results", {
    fit_m3 <- suppressMessages(suppressWarnings(
      .nlmixr(one.cmt, .datM3, est = "nlm", list(print = 0))
    ))
    fit_m4 <- suppressMessages(suppressWarnings(
      .nlmixr(one.cmt, .datM4, est = "nlm", list(print = 0))
    ))
    expect_false(isTRUE(all.equal(fit_m3$objf, fit_m4$objf)))
  })

})
