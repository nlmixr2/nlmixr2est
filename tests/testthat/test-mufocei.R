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
})
