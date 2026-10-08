nmTest({
  # TIME read in a PK-type statement: rxSolve(nonmem = TRUE) reads the record
  # time there, while d/dt() keeps the continuous time (#1167)
  .pkTimeModel <- function() {
    ini({
      tcl <- log(3)
      eta.cl ~ 0.1
      add.sd <- 0.1
    })
    model({
      cl <- exp(tcl + eta.cl) * (1 + (1 - exp(-0.05 * time)))
      d / dt(central) <- -cl / 30 * central
      cp <- central / 30
      cp ~ add(add.sd)
    })
  }

  .pkTimeData <- function() {
    # explicit doses: an ADDL dose expanded into $dataSav becomes a record
    .ev <- rxode2::et(amt = 100, time = c(0, 24, 48)) |>
      rxode2::et(c(1, 6, 12, 23, 30, 40, 47, 60, 72)) |>
      rxode2::et(id = 1:6)
    .p <- data.frame(id = 1:6, tcl = log(3), eta.cl = 0.1 * c(-2, -1, 0, 1, 2, 3), add.sd = 0)
    .s <- suppressWarnings(rxode2::rxSolve(.pkTimeModel, .ev, params = .p, nonmem = TRUE, covsInterpolation = "nocb", returnType = "data.frame"))
    .d <- as.data.frame(.ev)
    .d$dv <- NA_real_
    .d$dv[.d$evid == 0] <- .s$cp * (1 + 0.01 * sin(seq_along(.s$cp)))
    .d <- .d[, c("id", "time", "amt", "evid", "dv")]
    names(.d) <- toupper(names(.d))
    .d$AMT[is.na(.d$AMT)] <- 0
    .d
  }

  # IPRED and the fit's ETAs, re-solved directly by rxode2
  .pkTimeSolve <- function(fit, data, nonmem) {
    .p <- data.frame(id = fit$eta$ID, tcl = unname(fit$theta["tcl"]), eta.cl = fit$eta$eta.cl, add.sd = 0)
    suppressWarnings(rxode2::rxSolve(.pkTimeModel, data, params = .p, nonmem = nonmem, covsInterpolation = "nocb", returnType = "data.frame"))$cp
  }

  test_that("rx_time_pk~t is moved from ..lhs0 to the generated-model prologue", {
    .e <- new.env()
    .e$..lhs0 <- c(rx_time_pk = "rx_time_pk~t", a = "a=1")
    .e$..mtime <- "mtime(m)~2"
    .rxPkTimeToPrologue(.e)
    expect_equal(.e$..lhs0, c(a = "a=1"))
    expect_equal(.e$..mtime, c("rx_time_pk~t", "mtime(m)~2"))
    .e <- new.env()
    .e$..lhs0 <- c(a = "a=1")
    .rxPkTimeToPrologue(.e)
    expect_null(.e$..mtime)
  })

  .pkTimeChk <- function(model) {
    .mv <- rxode2::rxModelVars(model)
    expect_true(grepl("rx_time_pk~t", rxode2::rxNorm(.mv), fixed = TRUE))
    expect_false("rx_time_pk" %in% .mv$params)
  }

  test_that("generated models define rx_time_pk instead of reading it as a parameter (#1167)", {
    skip_if_not(.rxSHasPkTime(), "rxode2 without rxS(pkTime=)")
    .ui <- rxode2::rxode2(.pkTimeModel)
    .chk <- .pkTimeChk
    .m <- .ui$foceiModel
    .chk(.m$inner)
    .chk(.m$predOnly)
    # the inlined clearance reads the PK time, not the continuous one
    expect_true(grepl("d/dt(central)=", rxode2::rxNorm(.m$inner), fixed = TRUE))
    expect_false(grepl("exp(-0.05*t)", rxode2::rxNorm(.m$inner), fixed = TRUE))
    .chk(.ui$nlmeRxModel)
    .chk(.ui$saemModel)

    # a model that reads time only in d/dt() is unchanged
    .ode <- function() {
      ini({
        tcl <- log(3)
        eta.cl ~ 0.1
        add.sd <- 0.1
      })
      model({
        d / dt(central) <- -exp(tcl + eta.cl) * (2 - exp(-0.05 * time)) / 30 * central
        cp <- central / 30
        cp ~ add(add.sd)
      })
    }
    expect_false(grepl("rx_time_pk", rxode2::rxNorm(rxode2::rxode2(.ode)$foceiModel$inner), fixed = TRUE))
  })

  test_that("focei with rxControl(nonmem=TRUE) matches rxSolve(nonmem=TRUE) (#1167)", {
    skip_if_not(.rxSHasPkTime(), "rxode2 without rxS(pkTime=)")
    .d <- .pkTimeData()
    .fit <- .nlmixr(
      .pkTimeModel, .d, "focei",
      foceiControl(
        rxControl = rxode2::rxControl(covsInterpolation = "nocb", nonmem = TRUE),
        maxOuterIterations = 0L, covMethod = "", print = 0
      )
    )
    expect_true(.residCovsInterpolation(.fit)$nonmem)
    expect_true(grepl("rx_time_pk~t", rxode2::rxNorm(.fit$foceiModel$inner), fixed = TRUE))
    expect_equal(.fit$IPRED, .pkTimeSolve(.fit, .d, TRUE), tolerance = 1e-4)
    expect_gt(max(abs(.fit$IPRED / .pkTimeSolve(.fit, .d, FALSE) - 1)), 0.01)

    # augPred() and vpcSim() solve with the fit's nonmem setting too
    .seen <- logical(0)
    .orig <- rxode2::rxSolve
    local_mocked_bindings(
      rxSolve = function(...) {
        .seen <<- c(.seen, isTRUE(list(...)$nonmem))
        .orig(...)
      },
      .package = "rxode2"
    )
    expect_true(nrow(suppressWarnings(augPred(.fit))) > 0)
    expect_true(length(.seen) > 0 && all(.seen))
    .seen <- logical(0)
    expect_true(nrow(suppressWarnings(vpcSim(.fit, n = 2))) > 0)
    expect_true(length(.seen) > 0 && all(.seen))

    .fit0 <- .nlmixr(
      .pkTimeModel, .d, "focei",
      foceiControl(
        rxControl = rxode2::rxControl(covsInterpolation = "nocb"),
        maxOuterIterations = 0L, covMethod = "", print = 0
      )
    )
    expect_null(.residCovsInterpolation(.fit0)$nonmem)
    expect_equal(.fit0$IPRED, .pkTimeSolve(.fit0, .d, FALSE), tolerance = 1e-4)
    # the inner problem solved with the record time too, not just the table
    expect_false(isTRUE(all.equal(.fit$objf, .fit0$objf, tolerance = 1e-3)))
  })

  test_that("the analytic outer gradient reads the PK time under nonmem=TRUE (#1167)", {
    skip_if_not(.rxSHasPkTime(), "rxode2 without rxS(pkTime=)")
    .ctl <- function(...) {
      foceiControl(
        rxControl = rxode2::rxControl(covsInterpolation = "nocb", nonmem = TRUE),
        maxOuterIterations = 0L, covMethod = "", print = 0, ...
      )
    }
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.pkTimeModel))
    rxode2::rxAssignControlValue(.ui, "fast", TRUE)
    .pkTimeChk(.ui$foceiOuter$augMod)
    .d <- .pkTimeData()
    .fit <- .nlmixr(.pkTimeModel, .d, "focei", .ctl(fast = TRUE))
    expect_gt(.fit$env$nAnalyticGradDirect, 0)
    .g <- .foceiGradDirect(.fit)
    .ofv <- function(val) .nlmixr(rxode2::ini(.pkTimeModel, tcl = val), .d, "focei", .ctl())$objf
    .h <- 3e-3
    .fd <- (.ofv(log(3) + .h) - .ofv(log(3) - .h)) / (2 * .h)
    expect_equal(unname(.g["tcl"]), .fd, tolerance = 0.02)
  })

  test_that("nlminb estimates under rxControl(nonmem=TRUE) (#1167)", {
    skip_if_not(.rxSHasPkTime(), "rxode2 without rxS(pkTime=)")
    .pop <- function() {
      ini({
        tcl <- 1.3
        add.sd <- 0.1
      })
      model({
        cl <- exp(tcl) * (1 + (1 - exp(-0.05 * time)))
        d / dt(central) <- -cl / 30 * central
        cp <- central / 30
        cp ~ add(add.sd)
      })
    }
    .ev <- rxode2::et(amt = 100, time = c(0, 24, 48)) |>
      rxode2::et(c(1, 6, 12, 23, 30, 40, 47, 60, 72)) |>
      rxode2::et(id = 1:3)
    .s <- suppressWarnings(rxode2::rxSolve(.pop, .ev, params = c(tcl = log(3), add.sd = 0), nonmem = TRUE, covsInterpolation = "nocb", returnType = "data.frame"))
    .d <- as.data.frame(.ev)
    .d$dv <- NA_real_
    .d$dv[.d$evid == 0] <- .s$cp
    .d <- .d[, c("id", "time", "amt", "evid", "dv")]
    names(.d) <- toupper(names(.d))
    .d$AMT[is.na(.d$AMT)] <- 0
    .popUi <- rxode2::rxode2(.pop)
    .pkTimeChk(.popUi$nlmRxModel$predOnly)
    .pkTimeChk(.popUi$nlmSensModel$thetaGrad)
    .pkTimeChk(.popUi$nlmSensModel$predOnly)
    .fit <- .nlmixr(.pop, .d, "nlminb", nlminbControl(rxControl = rxode2::rxControl(covsInterpolation = "nocb", nonmem = TRUE)))
    expect_equal(unname(.fit$theta["tcl"]), log(3), tolerance = 1e-3)
  })
})
