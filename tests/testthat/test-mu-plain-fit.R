nmTest({
  # End-to-end fits for plain (covariate-free) mu-referenced profiling: the
  # mu-FOCEI family excludes plain mu-ref thetas from the outer optimizer via
  # intercept-only regression groups (see test-mu-plain.R for the UI-only
  # unit tests).

  .ocmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  theo_sd <- nlmixr2data::theo_sd

  fitFocei <- .getCachedFit(
    name = "mu-plain-focei",
    fitFn = function() .nlmixr(.ocmt, theo_sd, "focei", foceiControl(print = 0)),
    cacheFile = "fit-mu-plain-focei.rds"
  )

  test_that("irlsfocei profiles plain mu thetas out of the outer optimizer", {
    fit <- .getCachedFit(
      name = "mu-plain-irlsfocei",
      fitFn = function() .nlmixr(.ocmt, theo_sd, "irlsfocei",
                                 irlsfoceiControl(print = 0)),
      cacheFile = "fit-mu-plain-irlsfocei.rds"
    )
    # same optimum as plain focei (different path)
    expect_equal(unname(fit$theta), unname(fitFocei$theta), tolerance = 0.05)
    expect_equal(fit$objf, fitFocei$objf, tolerance = 0.5)
    # SEs come from the full-model covariance recompute
    .pf <- fit$parFixed
    for (.p in c("tka", "tcl", "tv")) {
      expect_false(is.na(suppressWarnings(as.numeric(.pf[.p, "SE"]))))
    }
    # only add.sd + the omegas were outer-optimized
    expect_equal(
      nlmixr2est:::.foceiMuSkipThetaNames(fit$ui,
        fit$ui$iniDf$name[!is.na(fit$ui$iniDf$ntheta)]),
      c("tka", "tcl", "tv"))
  })

  test_that("mufocei profiles plain mu thetas out of the outer optimizer", {
    fit <- .getCachedFit(
      name = "mu-plain-mufocei",
      fitFn = function() .nlmixr(.ocmt, theo_sd, "mufocei",
                                 mufoceiControl(print = 0)),
      cacheFile = "fit-mu-plain-mufocei.rds"
    )
    expect_equal(unname(fit$theta), unname(fitFocei$theta), tolerance = 0.05)
    expect_equal(fit$objf, fitFocei$objf, tolerance = 0.5)
  })

  test_that("irlsfoceif consumes the analytic gradient on the plain-profiled set", {
    fit <- .getCachedFit(
      name = "mu-plain-irlsfoceif",
      fitFn = function() .nlmixr(.ocmt, theo_sd, "irlsfoceif",
                                 irlsfoceiControl(print = 1)),
      cacheFile = "fit-mu-plain-irlsfoceif.rds"
    )
    expect_equal(unname(fit$theta), unname(fitFocei$theta), tolerance = 0.05)
    expect_equal(fit$objf, fitFocei$objf, tolerance = 0.5)
    expect_true("Analytic Gradient" %in% fit$parHistData$type)
  })

  test_that("iteration print shows plain mu thetas only in the mu rows", {
    out <- capture.output({
      nlmixr2est::nlmixr(.ocmt, theo_sd, "irlsfocei",
                         irlsfoceiControl(print = 1, maxOuterIterations = 2,
                                          covMethod = "", calcTables = FALSE))
    })
    muValueRows <- grep("^\\|   mu\\|.*tcl:\\s*[-0-9]", out, value = TRUE)
    expect_true(length(muValueRows) > 0)
    expect_true(all(grepl("tka:\\s*[-0-9]", muValueRows)))
    expect_true(all(grepl("tv:\\s*[-0-9]", muValueRows)))
    headerRows <- grep("^\\|    #\\|", out, value = TRUE)
    expect_true(length(headerRows) > 0)
    expect_false(any(grepl("\\btka\\b|\\btcl\\b|\\btv\\b", headerRows)))
    expect_true(any(grepl("add\\.sd", headerRows)))
  })
})
