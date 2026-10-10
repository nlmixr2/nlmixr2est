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
})
