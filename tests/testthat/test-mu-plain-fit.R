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

  .ocmtBnd <- function() {
    ini({
      tka <- 0.45
      tcl <- c(0, 1, 4.7)
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
  .ocmtClamp <- function() {
    ini({
      tka <- c(-2, 0.1, 0.2)
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

  test_that("a bounded mu theta is profiled (clamped regression), inactive bound matches unbounded", {
    fit <- .getCachedFit(
      name = "mu-plain-irls-bounded",
      fitFn = function() .nlmixr(.ocmtBnd, theo_sd, "irlsfocei",
                                 irlsfoceiControl(print = 0, covMethod = "",
                                                  calcTables = FALSE)),
      cacheFile = "fit-mu-plain-irls-bounded.rds"
    )
    fitFree <- .getCachedFit(
      name = "mu-plain-irls-free",
      fitFn = function() .nlmixr(.ocmt, theo_sd, "irlsfocei",
                                 irlsfoceiControl(print = 0, covMethod = "",
                                                  calcTables = FALSE)),
      cacheFile = "fit-mu-plain-irls-free.rds"
    )
    # tcl is profiled out despite its bounds
    expect_true("tcl" %in%
                  nlmixr2est:::.foceiMuSkipThetaNames(
                    fit$ui, fit$ui$iniDf$name[!is.na(fit$ui$iniDf$ntheta)]))
    # the interior optimum is unaffected by the inactive bound
    expect_equal(unname(fit$theta), unname(fitFree$theta), tolerance = 1e-4)
    expect_equal(fit$objf, fitFree$objf, tolerance = 1e-4)
  })

  test_that("an active bound clamps the regression update and is reported once", {
    fit <- .getCachedFit(
      name = "mu-plain-irls-clamped",
      fitFn = function() .nlmixr(.ocmtClamp, theo_sd, "irlsfocei",
                                 irlsfoceiControl(print = 0, covMethod = "",
                                                  calcTables = FALSE)),
      cacheFile = "fit-mu-plain-irls-clamped.rds"
    )
    # pinned exactly at the upper bound
    expect_equal(unname(fit$theta["tka"]), 0.2)
    # exactly ONE clamp note naming tka, marked at-bound
    .notes <- grep("mu-referenced regression clamped", fit$runInfo, value = TRUE)
    expect_length(.notes, 1L)
    expect_true(grepl("'tka'", .notes))
    expect_true(grepl("final estimate at bound", .notes))
    # plain focei with the same bound ends at the same place
    fitFocei2 <- .getCachedFit(
      name = "mu-plain-focei-clamped",
      fitFn = function() .nlmixr(.ocmtClamp, theo_sd, "focei",
                                 foceiControl(print = 0, covMethod = "",
                                              calcTables = FALSE)),
      cacheFile = "fit-mu-plain-focei-clamped.rds"
    )
    expect_equal(unname(fitFocei2$theta["tka"]), 0.2, tolerance = 1e-3)
    expect_equal(fit$objf, fitFocei2$objf, tolerance = 0.5)
    # a single-pass clamp cap still yields a feasible (in-bounds) fit
    fit1 <- .getCachedFit(
      name = "mu-plain-irls-clamp1",
      fitFn = function() .nlmixr(.ocmtClamp, theo_sd, "irlsfocei",
                                 irlsfoceiControl(print = 0, covMethod = "",
                                                  calcTables = FALSE,
                                                  muModelClampRetries = 1L)),
      cacheFile = "fit-mu-plain-irls-clamp1.rds"
    )
    expect_true(unname(fit1$theta["tka"]) <= 0.2)
    expect_true(unname(fit1$theta["tka"]) >= -2)
  })
})
