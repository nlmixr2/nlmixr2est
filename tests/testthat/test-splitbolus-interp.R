nmTest({
  # splitBolus() is applied once, when .foceiPreProcessData() translates the data
  # with the user's normModel.  The generated estimation models must therefore NOT
  # declare splitBolus() themselves, or the already-split doses get split twice.
  .mkSplit <- function(split) {
    .bdy <- c("ka <- exp(tka + eta.ka)", "cl <- exp(tcl)", "v <- exp(tv)",
              "d/dt(depot) <- -ka * depot",
              "d/dt(central) <- ka * depot - cl / v * central",
              if (split) "splitBolus(depot, depot, central)",
              "cp <- central / v", "cp ~ add(add.sd)")
    eval(parse(text=paste0(
      "function() {\n ini({tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.1; add.sd <- 0.7})\n",
      " model({\n", paste(.bdy, collapse="\n"), "\n })\n}")))
  }

  .d <- nlmixr2data::theo_sd
  .d <- .d[.d$ID <= 3, ]
  .d$EVID[.d$AMT > 0] <- 1L
  # equivalent of splitBolus(depot, depot, central): the same bolus into central
  .dose <- .d[.d$EVID == 1, ]
  .dose$CMT <- 2
  .dExpl <- rbind(.d, .dose)
  .dExpl <- .dExpl[order(.dExpl$ID, .dExpl$TIME, -.dExpl$EVID), ]

  test_that("generated models do not re-emit splitBolus()", {
    .ui <- .mkSplit(TRUE)()
    .ui <- rxode2::rxode2(.ui)
    expect_true(length(rxode2::rxModelVars(.ui)$splitBolus) > 0L)
    .noSplit <- function(mod) {
      expect_length(rxode2::rxModelVars(mod)$splitBolus, 0L)
    }
    .noSplit(rxUiGet.saemModel(list(.ui)))
    .noSplit(rxUiGet.saemModelPred(list(.ui))$predOnly)
    .noSplit(rxUiGet.nlmeRxModelFD(list(.ui)))
    .fm <- rxUiGet.foceiModel(list(.ui))
    for (.n in c("inner", "predOnly", "predNoLhs")) .noSplit(.fm[[.n]])
  })

  test_that("splitBolus() doses are split exactly once during estimation", {
    .ctl <- foceiControl(maxOuterIterations=0, maxInnerIterations=20, covMethod="",
                         print=0, calcTables=FALSE)
    .split <- nlmixr2(.mkSplit(TRUE), .d, est="focei", control=.ctl)
    .expl <- nlmixr2(.mkSplit(FALSE), .dExpl, est="focei", control=.ctl)
    .none <- nlmixr2(.mkSplit(FALSE), .d, est="focei", control=.ctl)
    expect_equal(.split$objf, .expl$objf, tolerance=1e-6)
    # guard against the split being dropped entirely
    expect_true(abs(.split$objf - .none$objf) > 1)
  })

  # est="nlm"/"nls" build their models from symengine parts, which drop the
  # covariate interpolation declared in the model
  .mkInterp <- function(interp) {
    .bdy <- c(if (!is.null(interp)) paste0(interp, "(WT)"),
              "ka <- exp(tka)", "cl <- exp(tcl + 0.05 * WT)", "v <- exp(tv)",
              "d/dt(depot) <- -ka * depot",
              "d/dt(central) <- ka * depot - cl / v * central",
              "cp <- central / v", "cp ~ add(add.sd)")
    eval(parse(text=paste0(
      "function() {\n ini({tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7})\n",
      " model({\n", paste(.bdy, collapse="\n"), "\n })\n}")))
  }

  test_that("nlm/nls honor the model's covariate interpolation", {
    .di <- .d
    withr::with_seed(1, .di$WT <- round(runif(nrow(.di), 60, 100), 1))
    .objf <- function(est, ctl, f) nlmixr2(f, .di, est=est, control=ctl)$objf
    .nctl <- nlmControl(print=0, calcTables=FALSE, iterlim=1)
    expect_false(isTRUE(all.equal(.objf("nlm", .nctl, .mkInterp("locf")),
                                  .objf("nlm", .nctl, .mkInterp("nocb")))))
    .sctl <- nlsControl(maxiter=1)
    expect_false(isTRUE(all.equal(.objf("nls", .sctl, .mkInterp("locf")),
                                  .objf("nls", .sctl, .mkInterp("nocb")))))
    # the interpolation reaches the compiled gradient model, not just predOnly
    .ui <- rxode2::rxode2(.mkInterp("nocb"))
    .interp <- rxode2::rxModelVars(rxUiGet.nlmSensModel(list(.ui))$thetaGrad)$interp
    expect_true(any(.interp != 0))
  })
})
