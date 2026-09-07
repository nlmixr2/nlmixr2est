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

  # Sets up (and tears down) an nlm environment for .mod with the given
  # normType, and returns the scaled starting parameter vector -- the same
  # quantity trust/nlm/bobyqa/etc. actually optimize over.
  .scaledParFor <- function(normType) {
    .ui <- rxode2::rxode2(.mod)
    .dat <- nlmixr2data::theo_sd
    .ctl <- nlmControl(print = 0, normType = normType, scaleType = "nlmixr2",
                       calcTables = FALSE, iterlim = 1)
    .ret <- new.env(parent = emptyenv())
    .foceiPreProcessData(.dat, .ret, .ui, .ctl$rxControl)
    .p <- setNames(.ui$nlmParIni, .ui$nlmParName)
    on.exit(.nlmFreeEnv())
    .nlmSetupEnv(.p, .ui, .ret$dataSav, .ui$nlmSensModel, .ctl)
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
})
