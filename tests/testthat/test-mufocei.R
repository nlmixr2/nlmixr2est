nmTest({
  test_that("mufocei recovers comparable estimates to plain focei", {
    # Local model/fits (not the shared helper-zzz-fits.R cache): this test
    # is specifically about the mu-referenced FOCEI restart-loop engine
    # itself, per the "WHEN TO KEEP LOCAL FITS" guidance in
    # helper-zzz-fits.R.

    theo_sd2 <- nlmixr2data::theo_sd
    theo_sd2$logWT <- log(theo_sd2$WT / 70)

    mod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        allo.cl <- 0.75 # mu-ref covariate coefficient (has an eta on cl)
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + allo.cl * logWT)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    fitFocei <- .getCachedFit(
      name = "mufocei-compare-focei",
      fitFn = function() .nlmixr(mod, theo_sd2, "focei", foceiControl(print = 0)),
      cacheFile = "fit-mufocei-compare-focei.rds"
    )

    fitMu <- .getCachedFit(
      name = "mufocei-basic",
      fitFn = function() .nlmixr(mod, theo_sd2, "mufocei", mufoceiControl(print = 0)),
      cacheFile = "fit-mufocei-basic.rds"
    )

    # structural: a normal, fully-classed FOCEI-family fit object
    expect_true(inherits(fitMu, "nlmixr2FitData"))
    expect_true(inherits(fitMu, "nlmixr2FitCore"))

    # the mu-ref-covariate theta/its coefficient get real standard errors
    # from the finalize step (proves the "unlock and refit" step ran)
    .pf <- fitMu$parFixed
    expect_true("allo.cl" %in% rownames(.pf))
    expect_false(is.na(suppressWarnings(as.numeric(.pf["allo.cl", "SE"]))))

    # recovered population theta / covariate coefficient / omega are
    # comparable to a plain focei fit of the same model on the same data
    .thFocei <- fitFocei$theta
    .thMu <- fitMu$theta
    expect_equal(unname(.thMu["tka"]), unname(.thFocei["tka"]), tolerance = 0.05)
    expect_equal(unname(.thMu["tcl"]), unname(.thFocei["tcl"]), tolerance = 0.05)
    expect_equal(unname(.thMu["tv"]), unname(.thFocei["tv"]), tolerance = 0.05)
    expect_equal(unname(.thMu["allo.cl"]), unname(.thFocei["allo.cl"]), tolerance = 0.2)
    expect_equal(unname(.thMu["add.sd"]), unname(.thFocei["add.sd"]), tolerance = 0.05)

    .omFocei <- diag(fitFocei$omega)
    .omMu <- diag(fitMu$omega)
    expect_equal(.omMu, .omFocei, tolerance = 0.2)

    # objective function values should be in the same ballpark (not an
    # exact match -- different optimization paths)
    expect_equal(fitMu$objf, fitFocei$objf, tolerance = 0.1)
  })

  test_that("mufocei respects a user-fixed covariate coefficient", {
    theo_sd2 <- nlmixr2data::theo_sd
    theo_sd2$logWT <- log(theo_sd2$WT / 70)

    modFixed <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        allo.cl <- fix(0.75)
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + allo.cl * logWT)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    fitFixed <- .getCachedFit(
      name = "mufocei-fixed-coef",
      fitFn = function() .nlmixr(modFixed, theo_sd2, "mufocei", mufoceiControl(print = 0)),
      cacheFile = "fit-mufocei-fixed-coef.rds"
    )

    expect_equal(unname(fitFixed$theta["allo.cl"]), 0.75, tolerance = 1e-6)
  })

  test_that("mufocei's live iteration print shows mu-group theta values, blank on gradient rows", {
    # Phase 5: the mu-group population/covariate thetas are excluded from
    # the outer optimizer's own parameter vector, so the shared
    # scale.h-based per-iteration table never shows them. printMuGroupThetaRow()
    # (src/inner.cpp) adds one extra "|   mu|..." row after each real
    # parameter print (current regression-updated value) and after each
    # gradient print (blank/NA, since these thetas are never part of the
    # gradient finite-difference -- there is no restart-loop involved,
    # the value comes from updateMuGroups() running inside innerOpt()).
    #
    # capture.output() must wrap the *un-suppressed* nlmixr() call directly
    # (not the .nlmixr()/suppressMessages() test helper) -- this console
    # output is otherwise swallowed by suppressMessages() in this R
    # session, unrelated to whether printMuGroupThetaRow() itself ran.
    theo_sd2 <- nlmixr2data::theo_sd
    theo_sd2$logWT <- log(theo_sd2$WT / 70)

    mod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        allo.cl <- 0.75
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + allo.cl * logWT)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    out <- capture.output({
      nlmixr2est::nlmixr(mod, theo_sd2, "mufocei",
                          mufoceiControl(print = 1, maxOuterIterations = 2))
    })

    muValueRows <- grep("^\\|   mu\\|.*tcl:\\s*[-0-9]", out, value = TRUE)
    muNaRows <- grep("^\\|   mu\\|.*tcl:\\s*NA", out, value = TRUE)
    gradRows <- grep("^\\|    [GFCMS]\\|", out, value = TRUE)

    expect_true(length(muValueRows) > 0)
    expect_true(length(muNaRows) > 0)
    # every mu-group value row also shows allo.cl (the covariate
    # coefficient), and every blank row shows NA for it too
    expect_true(all(grepl("allo\\.cl:\\s*[-0-9]", muValueRows)))
    expect_true(all(grepl("allo\\.cl:\\s*NA", muNaRows)))
    # a standard (non-mu-group) theta still shows a real numeric gradient
    # on the same gradient-print events
    expect_true(length(gradRows) > 0)
    expect_true(any(grepl("[-0-9]\\.[0-9]", gradRows)))
  })
})
