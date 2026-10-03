nmTest({
  # Issue #995: scaleSetup()'s per-normType setup loops in src/scale.h used
  # `for (unsigned int k = scale->npars-1; k--;)`, which evaluates the
  # PRE-decrement k for truthiness before the body runs, so the body never
  # executes with k=npars-1 -- the top-indexed parameter was silently
  # dropped from the mean/std/len accumulators (and its scaleC never reset
  # to NA_REAL). normType="mean"/"std"/"len" are reachable with the default
  # scaleType="nlmixr2" (nlmControl()'s own default), so this is a real,
  # user-visible numeric bug, not just an internal accounting detail.
  .mod <- function() {
    ini({
      tka <- 0.5
      tcl <- 1.2
      tv <- 3.0
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka); cl <- exp(tcl); v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  }

  # Loads an nlm problem for .mod under `ctl` and returns its starting
  # parameters; the caller frees it with .nlmFreeEnv().
  .loadMod <- function(ctl, model = .mod) {
    .ui <- rxode2::rxode2(model)
    .ret <- new.env(parent = emptyenv())
    .foceiPreProcessData(nlmixr2data::theo_sd, .ret, .ui, ctl$rxControl)
    .p <- setNames(.ui$nlmParIni, .ui$nlmParName)
    .nlmSetupEnv(.p, .ui, .ret$dataSav, .ui$nlmSensModel, ctl)
    .p
  }

  # The scaled starting parameter vector for the given normType -- the same
  # quantity trust/nlm/bobyqa/etc. actually optimize over.
  .scaledParFor <- function(normType) {
    on.exit(.nlmFreeEnv())
    .p <- .loadMod(nlmControl(print = 0, normType = normType, scaleType = "nlmixr2", calcTables = FALSE, iterlim = 1))
    list(par = .p, scaled = nlmScalePar(.p))
  }

  test_that("normType='mean' includes every parameter, not just the first n-1 (#995)", {
    skip_on_cran()
    .r <- .scaledParFor("mean")
    .want <- (.r$par - mean(.r$par)) / (max(.r$par) - min(.r$par))
    expect_equal(unname(.r$scaled), unname(.want), tolerance = 1e-8)
  })

  test_that("normType='std' includes every parameter, not just the first n-1 (#995)", {
    skip_on_cran()
    .r <- .scaledParFor("std")
    .want <- (.r$par - mean(.r$par)) / sd(.r$par)
    expect_equal(unname(.r$scaled), unname(.want), tolerance = 1e-8)
  })

  test_that("normType='len' includes every parameter, not just the first n-1 (#995)", {
    skip_on_cran()
    .r <- .scaledParFor("len")
    .want <- .r$par / sqrt(sum(.r$par^2))
    expect_equal(unname(.r$scaled), unname(.want), tolerance = 1e-8)
  })

  # every theta at the same value v, scaled by scaleType = "norm": x = (p - c1)/c2
  .unscaledOnes <- function(v) {
    on.exit(.nlmFreeEnv())
    .m <- rxode2::ini(.mod, tka = v, tcl = v, tv = v, add.sd = v)
    .ctl <- nlmControl(print = 0, scaleType = "norm", calcTables = FALSE, iterlim = 1)
    .acc <- new.env(parent = emptyenv())
    .acc$warn <- character(0)
    .p <- withCallingHandlers(.loadMod(.ctl, .m), warning = function(w) {
      .acc$warn <- c(.acc$warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    list(p = .p, ones = nlmUnscalePar(rep(1, 4)), warn = .acc$warn)
  }

  test_that("equal initial estimates are normalized to unit length", {
    skip_on_cran()
    .r <- .unscaledOnes(0.7)
    expect_identical(.r$warn, "all parameters are the same value, switch to length normType")
    # c1 = 0, c2 = sqrt(4 * 0.49) = 1.4
    expect_equal(unname(.r$ones), rep(1.4, 4), tolerance = 1e-12)
  })

  test_that("all-zero initial estimates are not normalized", {
    skip_on_cran()
    .r <- .unscaledOnes(0)
    expect_identical(.r$warn, "all parameters are zero, cannot scale, run unscaled")
    # c1 = 0, c2 = 1
    expect_equal(unname(.r$ones), rep(1, 4))
  })

  test_that("scaleType='mult' scales and unscales only for a positive scaleTo", {
    skip_on_cran()
    # scaleTo <= 0 means no scaling; nlmControl() rejects a negative scaleTo,
    # so it is set on the built control
    .ctl <- nlmControl(print = 0, scaleType = "mult", calcTables = FALSE, iterlim = 1)
    .ctl$scaleTo <- -1
    on.exit(.nlmFreeEnv())
    .p <- .loadMod(.ctl)
    expect_identical(nlmScalePar(.p), unname(.p))
    expect_identical(nlmUnscalePar(.p), .p)
    expect_identical(.nlmAdjustCov(diag(4), .p), diag(4))
  })

  # FOCEi's outer problem normalizes its parameters (the thetas, then the omega
  # parameters) by the same rule: the mean, sd or length covers every one of
  # them, the last omega parameter included
  .modEta <- function() {
    ini({
      tka <- 0.5
      tcl <- 1.2
      tv <- 3.0
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  }

  # FOCEi's starting parameters and their scaled values, from the first
  # evaluation of the fit's parameter history
  .foceiScaledFor <- function(normType) {
    .f <- suppressMessages(suppressWarnings(nlmixr(
      .modEta,
      nlmixr2data::theo_sd,
      "focei",
      control = foceiControl(
        print = 0,
        maxOuterIterations = 1L,
        covMethod = "",
        calcTables = FALSE,
        outerOpt = "lbfgsb3c",
        normType = normType
      )
    )))
    .ph <- .f$parHistData
    .p <- setdiff(names(.ph), c("iter", "type", "objf"))
    list(
      par = unlist(.ph[.ph$type == "Unscaled", .p][1, ]),
      scaled = unlist(.ph[.ph$type == "Scaled", .p][1, ])
    )
  }

  # scaled = (par - c1) / c2 for each normalization
  .normWant <- list(
    rescale2 = function(p) (p - (max(p) + min(p)) / 2) / ((max(p) - min(p)) / 2),
    rescale = function(p) (p - min(p)) / (max(p) - min(p)),
    mean = function(p) (p - mean(p)) / (max(p) - min(p)),
    std = function(p) (p - mean(p)) / sd(p),
    len = function(p) p / sqrt(sum(p^2))
  )

  test_that("FOCEi normalizes over every parameter for each normType", {
    skip_on_cran()
    for (.nt in names(.normWant)) {
      .r <- .foceiScaledFor(.nt)
      expect_equal(unname(.r$scaled), unname(.normWant[[.nt]](.r$par)), tolerance = 1e-8, label = .nt)
    }
  })
})
