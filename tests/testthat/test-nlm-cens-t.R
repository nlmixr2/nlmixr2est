nmTest({
  # Tests for M2/M3/M4 censoring support for t()/cauchy() endpoints on the
  # NLM estimation engine (#979).  Cauchy is Student-t with nu=1, so
  # doCensT1() serves both; verify that identity too.

  one.cmt.t <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      nu <- fix(5)
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) + t(nu)
    })
  }

  one.cmt.cauchy <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) + cauchy()
    })
  }

  one.cmt.t1 <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      nu <- fix(1)
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) + t(nu)
    })
  }

  .dat <- nlmixr2data::theo_sd
  .dat <- .dat[.dat$ID <= 4, ]

  .LLOQ <- 3.0
  .datM3 <- .dat
  .datM3$CENS <- ifelse(.datM3$DV < .LLOQ & .datM3$EVID == 0, 1L, 0L)
  .datM3$DV[.datM3$CENS == 1] <- .LLOQ

  # A "naive" comparison dataset: same DV values, but CENS all zero, so the
  # censored rows are scored as if they were real, uncensored observations
  # AT the LLOQ.  If the M2/M3/M4 correction is being applied, the fitted
  # objective must differ from this baseline.
  .datNaive <- .datM3
  .datNaive$CENS <- 0L

  test_that("nlm t() endpoint honors M3 censoring (#979)", {
    fitCens <- suppressMessages(suppressWarnings(
      .nlmixr(one.cmt.t, .datM3, est = "nlm", list(print = 0))))
    expect_equal(as.character(fitCens$censInformation), "M3 censoring")

    fitNaive <- suppressMessages(suppressWarnings(
      .nlmixr(one.cmt.t, .datNaive, est = "nlm", list(print = 0))))
    expect_equal(as.character(fitNaive$censInformation), "No censoring")

    # the M2/M3/M4 correction must change the objective relative to the
    # naive (censoring-ignored) baseline -- this is the exact silent-
    # corruption issue #979 reported.
    expect_true(abs(fitCens$objf - fitNaive$objf) > 1e-6)
  })

  test_that("nlm cauchy() endpoint honors M3 censoring (#979)", {
    fitCens <- suppressMessages(suppressWarnings(
      .nlmixr(one.cmt.cauchy, .datM3, est = "nlm", list(print = 0))))
    expect_equal(as.character(fitCens$censInformation), "M3 censoring")

    fitNaive <- suppressMessages(suppressWarnings(
      .nlmixr(one.cmt.cauchy, .datNaive, est = "nlm", list(print = 0))))
    expect_equal(as.character(fitNaive$censInformation), "No censoring")

    expect_true(abs(fitCens$objf - fitNaive$objf) > 1e-6)
  })

  test_that("t(nu=1) censored objective matches cauchy() censored objective (#979)", {
    # cauchy() is Student-t with nu=1 -- doCensT1() is a single function
    # serving both, so the fitted -2LL at the same starting values/data
    # should agree to solver tolerance.
    fitT1 <- suppressMessages(suppressWarnings(
      .nlmixr(one.cmt.t1, .datM3, est = "nlm", list(print = 0))))
    fitCauchy <- suppressMessages(suppressWarnings(
      .nlmixr(one.cmt.cauchy, .datM3, est = "nlm", list(print = 0))))
    expect_equal(fitT1$objf, fitCauchy$objf, tolerance = 1e-3)
  })

  test_that("nlm t()/cauchy() censoring does not warn (#979 catch-all)", {
    expect_no_warning(suppressMessages(
      .nlmixr(one.cmt.t, .datM3, est = "nlm", list(print = 0))))
    expect_no_warning(suppressMessages(
      .nlmixr(one.cmt.cauchy, .datM3, est = "nlm", list(print = 0))))
  })
})
