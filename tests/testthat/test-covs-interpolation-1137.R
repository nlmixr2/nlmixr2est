nmTest({
  # A time-varying covariate, so locf and nocb give different predictions
  .covsIntModel <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      cl.crcl <- 0.5
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl) * (CRCL / 100)^cl.crcl
      v <- exp(tv + eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  }
  .covsIntData <- nlmixr2data::theo_sd
  .covsIntData$CRCL <- 90 + 2 * .covsIntData$TIME + .covsIntData$ID

  # Solve the fit's model at its thetas/etas; return cp matched to the table rows
  .covsIntSolve <- function(fit, covsInterpolation) {
    .p <- cbind(fit$eta, as.data.frame(as.list(fit$theta)))
    .ev <- .covsIntData
    names(.ev) <- tolower(names(.ev))
    names(.ev)[names(.ev) == "crcl"] <- "CRCL"
    .s <- suppressWarnings(rxode2::rxSolve(
      fit$ui$simulationModel,
      params = .p,
      events = .ev,
      covsInterpolation = covsInterpolation,
      returnType = "data.frame"
    ))
    .fd <- as.data.frame(fit)
    .m <- merge(
      data.frame(ID = as.integer(.fd$ID), TIME = .fd$TIME, IPRED = .fd$IPRED),
      data.frame(ID = as.integer(.s$id), TIME = .s$time, cp = .s$cp)
    )
    .m[.m$cp > 0, ]
  }

  test_that("the fit's table uses rxControl(covsInterpolation=) (#1137)", {
    for (.ci in c("locf", "nocb", "midpoint")) {
      .fit <- .nlmixr(
        .covsIntModel,
        .covsIntData,
        "focei",
        foceiControl(
          rxControl = rxode2::rxControl(covsInterpolation = .ci),
          maxOuterIterations = 0L,
          print = 0
        )
      )
      expect_equal(
        .residCovsInterpolation(.fit)$covsInterpolation,
        rxode2::rxControl(covsInterpolation = .ci)$covsInterpolation,
        ignore_attr = TRUE
      )
      .other <- if (.ci == "locf") "nocb" else "locf"
      .m <- .covsIntSolve(.fit, .ci)
      .mOther <- .covsIntSolve(.fit, .other)
      # matches the requested interpolation, not the other one
      expect_equal(.m$IPRED, .m$cp, tolerance = 1e-4)
      expect_gt(median(abs(.mOther$IPRED - .mOther$cp) / .mOther$cp), 1e-3)

      # augPred() defaults to the fit's interpolation too
      .ap <- augPred(.fit)
      .ap <- .ap[.ap$ind == "Individual", ]
      .ma <- merge(
        data.frame(id = as.integer(.m$ID), time = .m$TIME, IPRED = .m$IPRED),
        data.frame(id = as.integer(.ap$id), time = .ap$time, ap = .ap$values)
      )
      expect_true(nrow(.ma) > 0)
      # augPred's added grid rows move the midpoints, so midpoint is close only
      expect_equal(.ma$ap, .ma$IPRED, tolerance = if (.ci == "midpoint") 1e-3 else 1e-4)

      # vpcSim() (and so npde) simulates with it unless overridden
      .vDefault <- vpcSim(.fit, n = 2, seed = 42)
      .vSame <- vpcSim(.fit, n = 2, seed = 42, covsInterpolation = .ci)
      .vOther <- vpcSim(.fit, n = 2, seed = 42, covsInterpolation = .other)
      expect_equal(.vDefault$ipred, .vSame$ipred)
      expect_false(isTRUE(all.equal(.vDefault$ipred, .vOther$ipred)))
      # a setting simInfo already holds is replaced, not passed twice
      .vEvents <- vpcSim(.fit, n = 2, seed = 42, events = .fit$origData)
      expect_equal(.vEvents$ipred, .vDefault$ipred)
    }
  })
  test_that("the fit's table uses rxControl(naInterpolation=) (#1137)", {
    # naInterpolation fills NA covariates under linear/midpoint interpolation;
    # the table's covariate column shows which side it filled from
    .d <- .covsIntData
    .d$CRCL[.d$TIME > 0 & .d$TIME < 10] <- NA
    .naFit <- function(naInterpolation) {
      .nlmixr(
        .covsIntModel,
        .d,
        "focei",
        foceiControl(
          rxControl = rxode2::rxControl(
            covsInterpolation = "linear",
            naInterpolation = naInterpolation
          ),
          maxOuterIterations = 0L,
          print = 0
        )
      )
    }
    # CRCL of the first subject's rows inside the NA window
    .window <- function(df, id, time) {
      df$CRCL[id == id[1] & time > 0 & time < 10]
    }
    .fitNocb <- .naFit("nocb")
    .fdNocb <- as.data.frame(.fitNocb)
    .fdLocf <- as.data.frame(.naFit("locf"))
    .nocb <- .window(.fdNocb, .fdNocb$ID, .fdNocb$TIME)
    .locf <- .window(.fdLocf, .fdLocf$ID, .fdLocf$TIME)
    .first <- .covsIntData[.covsIntData$ID == .covsIntData$ID[1], ]
    .before <- .first$CRCL[.first$TIME == 0][1]
    .after <- .first$CRCL[.first$TIME >= 10][1]
    expect_true(length(.nocb) > 0)
    expect_true(all(.locf == .before))
    expect_true(all(.nocb == .after))

    # vpcSim() fills them the same way unless overridden
    .v <- vpcSim(.fitNocb, n = 1, seed = 42)
    .vLocf <- vpcSim(.fitNocb, n = 1, seed = 42, naInterpolation = "locf")
    expect_true(all(.window(.v, .v$id, .v$time) == .after))
    expect_true(all(.window(.vLocf, .vLocf$id, .vLocf$time) == .before))
  })
})
