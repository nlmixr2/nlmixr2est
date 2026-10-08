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
  .loadMod <- function(ctl) {
    .ui <- rxode2::rxode2(.mod)
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

  test_that("a derivative-based scaleC is band-guarded in R and in C++ (#994)", {
    skip_on_cran()
    # tcl's near-zero starting gradient sends |gradTo/gradient| far out of band
    .bounded <- function() {
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
    .ui <- rxode2::rxode2(.bounded)
    .ctl <- nlmControl(print = 0, scaleType = "nlmixr2", calcTables = FALSE, iterlim = 1)
    .ret <- new.env(parent = emptyenv())
    .foceiPreProcessData(nlmixr2data::theo_sd, .ret, .ui, .ctl$rxControl)
    .p <- setNames(.ui$nlmParIni, .ui$nlmParName)
    on.exit(.nlmFreeEnv())
    .env <- .nlmSetupEnv(.p, .ui, .ret$dataSav, .ui$nlmSensModel, .ctl)
    # with scaleType="nlmixr2", du/dx is the scaleC the C++ side holds
    .used <- unname(nlmUnscalePar(.env$par.ini + 1) - .p)
    .raw <- nlmGetScaleC(.p, .ctl$gradTo)
    expect_true(any(.raw < 0.1 | .raw > 10))
    .want <- mapply(.guardScaleC, .raw, .ui$scaleCtheta, USE.NAMES = FALSE)
    expect_equal(.env$scaleC, .want)
    expect_equal(.used, .want, tolerance = 1e-8)
  })
})
