# ar() correlations that are not an "ar" row of the fit's $iniDf: a fixed
# one and a modeled one
.arModFixedCor <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    eta.cl ~ 0.1
    add.sd <- 0.7
    ar1.cor <- fix(0.9)
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd) + ar(ar1.cor)
  })
}

.arModModeledCor <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    eta.cl ~ 0.1
    add.sd <- 0.7
    tcor <- 2
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    corv <- expit(tcor)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd) + ar(corv)
  })
}

# the ui as a focei fit builds from it, after literalFixRes writes a fixed
# residual parameter's value into the model
.arUiAsFit <- function(mod) {
  .ui <- if (is.function(mod)) rxode2::rxode2(mod) else mod
  .env <- nlmixr2global$nlmixr2EstEnv
  .old <- .env$uiUnfix
  on.exit(assign("uiUnfix", .old, envir = .env))
  .nlmixrPreprocessLiteralFix(.ui, "focei", NULL, foceiControl())$ui
}

nmTest({
  # ar() residuals are out of analytic outer-gradient scope (#1195)
  .arMod <- function() {
    ini({
      tka <- log(1.2)
      tcl <- log(0.2)
      tv <- log(5)
      eta.cl ~ 0.09
      add.sd <- 0.4
      ar1.cor <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd) + ar(ar1.cor)
    })
  }

  test_that("ar() residuals are detected and declined by the analytic gradient gates", {
    .ui <- rxode2::rxode2(.arMod)
    expect_true(.foceiUsesAr(.ui))
    expect_false(.foceiUsesAr(rxode2::rxode2(.arMod) |> rxode2::model(cp ~ add(add.sd))))
    expect_null(.foceiOuterDirs(.ui, "focei"))
    expect_false(.foceiLLGradInScope(.ui, "focei"))
  })

  test_that(".foceiUsesAr() sees every ar() residual on the ui the fit builds from", {
    .none <- rxode2::rxode2(.arMod) |> rxode2::model(cp ~ add(add.sd))
    expect_false(.foceiUsesAr(.arUiAsFit(.none)))
    expect_true(.foceiUsesAr(.arUiAsFit(.arMod)))
    # a fixed or modeled correlation is no "ar" row of $iniDf
    .fixed <- .arUiAsFit(.arModFixedCor)
    expect_false("ar" %in% .fixed$iniDf$err)
    expect_true(.foceiUsesAr(.fixed))
    .modeled <- .arUiAsFit(.arModModeledCor)
    expect_false("ar" %in% .modeled$iniDf$err)
    expect_true(.foceiUsesAr(.modeled))
  })

  test_that("fast = TRUE and the analytic covariance fall back to finite differences for ar()", {
    skip_on_cran()
    .ev <- rxode2::et(rxode2::et(amt = 100), seq(0.5, 24, by = 2))
    .sim <- rxode2::rxWithSeed(
      2026,
      rxode2::rxSolve(.arMod, .ev, nSub = 12, returnType = "data.frame", addDosing = TRUE)
    )
    .idc <- if ("sim.id" %in% names(.sim)) "sim.id" else "id"
    .dat <- data.frame(
      ID = .sim[[.idc]],
      TIME = .sim$time,
      EVID = ifelse(is.na(.sim$evid), 0, .sim$evid),
      AMT = ifelse(is.na(.sim$amt), 0, .sim$amt),
      DV = ifelse(is.na(.sim$evid) | .sim$evid == 0, .sim$sim, NA)
    )
    .dat$EVID[.dat$AMT > 0] <- 1
    .acc <- new.env(parent = emptyenv())
    .acc$msg <- character(0)
    .fit <- withCallingHandlers(
      suppressWarnings(nlmixr2(
        .arMod,
        .dat,
        "focei",
        foceiControl(fast = TRUE, print = 0L, maxOuterIterations = 2L, covMethod = "analytic", calcTables = FALSE)
      )),
      message = function(m) {
        .acc$msg <- c(.acc$msg, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
    expect_true(is.finite(.fit$objf))
    expect_true(any(grepl("ar() residual: the analytic 'fast' gradient does not apply", .acc$msg, fixed = TRUE)))
    expect_true(any(grepl("an ar() residual is out of analytic-covariance scope", .acc$msg, fixed = TRUE)))
    expect_false(.covBaseName(.fit$covMethod) == "analytic")
    expect_false(isTRUE(.fit$foceiControl$fast))
    expect_false(isTRUE(.fit$env$nAnalyticGradDirect > 0))
  })

  # the fixed correlation is no "ar" row of the ui FOCEi builds from; the
  # fast gradient used to be built anyway, without the AR(1) lines, and the
  # fit stopped ("required for solving: rx_arNf_cp, rx_arEp_cp, rx_arT_cp")
  test_that("fast = TRUE falls back for an ar() with a fixed or a modeled correlation", {
    skip_on_cran()
    for (.m in list(fixed = .arModFixedCor, modeled = .arModModeledCor)) {
      .msg <- testthat::capture_messages(
        .fit <- suppressWarnings(nlmixr2(
          .m,
          nlmixr2data::theo_sd,
          "focei",
          foceiControl(fast = TRUE, print = 0L, maxOuterIterations = 0L, covMethod = "", calcTables = FALSE)
        ))
      )
      expect_true(is.finite(.fit$objf))
      expect_true(any(grepl("ar() residual: the analytic 'fast' gradient does not apply", .msg, fixed = TRUE)))
      expect_false(isTRUE(.fit$foceiControl$fast))
    }
  })
})
