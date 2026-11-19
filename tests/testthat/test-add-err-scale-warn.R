nmTest({
  # #1383: warn when a combined error model's additive SD starts negligibly
  # small relative to the observed DV scale
  .mod <- function(addSd = 0.01, err = "comb", fixAdd = FALSE) {
    .f <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- 0.01
        prop.sd <- 0.2
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        cp <- linCmt()
        cp ~ add(add.sd) + prop(prop.sd)
      })
    }
    .ui <- rxode2::rxode2(.f)
    .ui <- rxode2::ini(.ui, add.sd = addSd)
    if (fixAdd) {
      .ui <- rxode2::ini(.ui, add.sd = fix(addSd))
    }
    if (err == "add") {
      .ui <- rxode2::model(.ui, cp ~ add(add.sd))
    } else if (err == "prop") {
      .ui <- rxode2::model(.ui, cp ~ prop(prop.sd))
    }
    .ui
  }

  .dat <- nlmixr2data::theo_sd
  .dat <- .dat[.dat$ID <= 4, ]
  .big <- .dat
  .big$DV <- .big$DV * 50 # median |DV| in the hundreds

  .warns <- function(ui, data, est = "focei") {
    .w <- character(0)
    withCallingHandlers(
      .preProcessAddErrScaleWarn(ui, est, data, list()),
      warning = function(w) {
        .w <<- c(.w, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    .w
  }

  test_that("tiny additive start on a combined model warns", {
    .w <- .warns(.mod(0.01), .big)
    expect_length(.w, 1L)
    expect_match(.w, "add sd 'add.sd'=0.01 tiny vs median |DV|", fixed = TRUE)
    expect_lt(nchar(.w), 75)
  })

  test_that("a data-scaled start does not warn", {
    .med <- stats::median(abs(.big$DV[.big$EVID == 0]))
    expect_length(.warns(.mod(0.1 * .med), .big), 0L)
  })

  test_that("additive-only, proportional-only and fixed additive do not warn", {
    expect_length(.warns(.mod(0.01, "add"), .big), 0L)
    expect_length(.warns(.mod(0.01, "prop"), .big), 0L)
    expect_length(.warns(.mod(0.01, fixAdd = TRUE), .big), 0L)
  })

  test_that("no positive DV and simulation do not warn", {
    .neg <- .big
    .neg$DV <- -abs(.neg$DV)
    expect_length(.warns(.mod(0.01), .neg), 0L)
    expect_length(.warns(.mod(0.01), .big, "rxSolve"), 0L)
  })

  test_that("multiple endpoints are checked per endpoint", {
    .f <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- 0.5
        prop.sd <- 0.2
        add.sd2 <- 0.01
        prop.sd2 <- 0.2
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        cp <- linCmt()
        cp2 <- cp
        cp ~ add(add.sd) + prop(prop.sd)
        cp2 ~ add(add.sd2) + prop(prop.sd2)
      })
    }
    .d1 <- .dat[.dat$EVID == 0, ]
    .d1$DVID <- 1L
    .d2 <- .d1
    .d2$DVID <- 2L
    .d2$DV <- .d2$DV * 50
    .dose <- .dat[.dat$EVID != 0, ]
    .dose$DVID <- NA_integer_
    .d <- rbind(.dose, .d1, .d2)
    .d <- .d[order(.d$ID, .d$TIME, -.d$EVID), ]
    .w <- .warns(rxode2::rxode2(.f), .d)
    expect_length(.w, 1L)
    expect_match(.w, "'add.sd2'", fixed = TRUE)
  })

  test_that("the warning reaches $runInfo and leaves the initial estimates alone", {
    skip_on_cran()
    .ui <- .mod(0.01)
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(.ui, .big, est = "focei", control = foceiControl(print = 0, maxOuterIterations = 0, covMethod = ""))
    ))
    expect_true(any(grepl("add sd 'add.sd'=0.01 tiny", fit$runInfo, fixed = TRUE)))
    # no outer iterations: the estimate is the (untouched) initial value
    expect_equal(unname(fit$theta["add.sd"]), 0.01)
  })
})
