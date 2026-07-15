nmTest({
  test_that("mu/irls covariance is the full base model at the converged point (ODE/matExp)", {
    # The mu/irls covariance step bails and is recomputed on the full base
    # focei/foce/focep model (muModel="none") AT the mu fit's converged theta + eta with
    # the inner problem FROZEN (maxOuter=maxInner=0) -- true to that point, not re-optimized.
    # It must equal the base method computed the SAME way at that point, and (analytic
    # default, well-identified) equal a normal base fit.
    theo_sd2 <- nlmixr2data::theo_sd

    odeNC <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; eta.cl ~ 0.3; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot; d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v; cp ~ add(add.sd) })
    }
    matNC <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; eta.cl ~ 0.3; add.sd <- 0.7 })
      model({ matExp()
        k_depot_central <- exp(tka + eta.ka)
        k_central_output <- exp(tcl + eta.cl) / exp(tv)
        cp <- central / exp(tv); cp ~ add(add.sd) })
    }

    # base method computed with the SAME frozen algorithm at a given theta/eta point
    .baseFrozen <- function(mfun, baseEst, baseCtl, ui, eta, cm) {
      .em <- as.matrix(eta[, setdiff(names(eta), "ID"), drop = FALSE])
      .nlmixr(rxode2::rxUiDecompress(ui), theo_sd2, baseEst,
              baseCtl(print = 0, covMethod = cm, maxOuterIterations = 0L,
                      maxInnerIterations = 0L, etaMat = .em))
    }

    .chk <- function(mfun, muEst, baseEst, muCtl, baseCtl) {
      # analytic is the default; the mu recompute freezes the inner problem, and for the
      # analytic observed information that reproduces the same value as a normal fit.
      .fM <- .nlmixr(mfun, theo_sd2, muEst, muCtl(print = 0))
      expect_equal(.fM$covMethod, "analytic")
      expect_false(is.na(suppressWarnings(as.numeric(.fM$parFixed["tcl", "SE"]))))
      # (a) equals the base method computed the SAME (frozen) way at the mu point
      .fF <- .baseFrozen(mfun, baseEst, baseCtl, .fM$ui, .fM$eta, "analytic")
      .sM <- sqrt(diag(.fM$cov)); .sF <- sqrt(diag(.fF$cov))
      .cmn <- intersect(names(.sM), names(.sF))
      expect_equal(unname(.sM[.cmn]), unname(.sF[.cmn]), tolerance = 1e-2)
      # (b) well-identified: also equals a normal base fit (converges to the same point)
      .fN <- .nlmixr(mfun, theo_sd2, baseEst, baseCtl(print = 0))
      .sN <- sqrt(diag(.fN$cov))
      expect_equal(unname(.sM[.cmn]), unname(.sN[.cmn]), tolerance = 5e-2)
    }

    for (.mod in list(odeNC, matNC)) {
      .chk(.mod, "mfocei",   "focei", mfoceiControl,   foceiControl)
      .chk(.mod, "ifocei", "focei", ifoceiControl, foceiControl)
      .chk(.mod, "mfoce",    "foce",  mfoceControl,    foceControl)
      .chk(.mod, "ifoce",  "foce",  ifoceControl,  foceControl)
      .chk(.mod, "mfocep",   "focep", mfocepControl,   focepControl)
      .chk(.mod, "ifocep", "focep", ifocepControl, focepControl)
    }
  })

  test_that("mu/irls covariance on a covariate (linear-component) model has finite SEs", {
    # The mu-ref/covariate ("linear") parameters (the ~8.5x-wrong case under the mu->phi
    # reduction) get finite SEs from the full-model recompute; the weakly-identified
    # covariate coefficient converges to a slightly different point than a normal fit, so
    # only its finiteness + the well-identified structural thetas are checked.
    theo_sd2 <- nlmixr2data::theo_sd
    theo_sd2$logWT <- log(theo_sd2$WT / 70)

    modOde <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; allo.cl <- 0.75
            eta.ka ~ 0.6; eta.cl ~ 0.3; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl + allo.cl * logWT); v <- exp(tv)
        d/dt(depot) <- -ka * depot; d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v; cp ~ add(add.sd) })
    }
    for (.est in c("mfocei", "ifocei")) {
      .ctl <- if (.est == "mfocei") mfoceiControl else ifoceiControl
      .fB <- .nlmixr(modOde, theo_sd2, "focei", foceiControl(print = 0))
      .fM <- .nlmixr(modOde, theo_sd2, .est, .ctl(print = 0))
      expect_equal(.fM$covMethod, "analytic")
      expect_false(is.na(suppressWarnings(as.numeric(.fM$parFixed["allo.cl", "SE"]))))
      expect_false(is.na(suppressWarnings(as.numeric(.fM$parFixed["tcl", "SE"]))))
      .sM <- sqrt(diag(.fM$cov)); .sB <- sqrt(diag(.fB$cov))
      .th <- c("tka", "tcl", "tv", "add.sd")
      expect_equal(unname(.sM[.th]), unname(.sB[.th]), tolerance = 0.1)
    }
  })

  test_that("mfocei recovers comparable estimates to plain focei", {
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
      name = "mfocei-compare-focei",
      fitFn = function() .nlmixr(mod, theo_sd2, "focei", foceiControl(print = 0)),
      cacheFile = "fit-mfocei-compare-focei.rds"
    )

    fitMu <- .getCachedFit(
      name = "mfocei-basic",
      fitFn = function() .nlmixr(mod, theo_sd2, "mfocei", mfoceiControl(print = 0)),
      cacheFile = "fit-mfocei-basic.rds"
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

    # mu thetas are recorded in the parameter history like plain focei
    expect_identical(names(fitMu$parHist), names(fitFocei$parHist))
    .u <- fitMu$parHistData[fitMu$parHistData$type == "Unscaled", ]
    expect_true(nrow(.u) > 0)
    expect_true(all(is.finite(.u$tcl)))
    expect_true(all(is.finite(.u$allo.cl)))
    expect_equal(unname(.u$tcl[nrow(.u)]), unname(.thMu["tcl"]), tolerance = 0.05)
    expect_equal(unname(.u$allo.cl[nrow(.u)]), unname(.thMu["allo.cl"]),
                 tolerance = 0.05)
  })

  test_that("mfocei respects a user-fixed covariate coefficient", {
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
      name = "mfocei-fixed-coef",
      fitFn = function() .nlmixr(modFixed, theo_sd2, "mfocei", mfoceiControl(print = 0)),
      cacheFile = "fit-mfocei-fixed-coef.rds"
    )

    expect_equal(unname(fitFixed$theta["allo.cl"]), 0.75, tolerance = 1e-6)
  })

  test_that("mfocei regression-updates a bounded mu-ref covariate coefficient with clamping", {
    # A bounded covariate coefficient is regression-updated like any other
    # (the update is clamped to the bounds by the active-set loop in
    # updateMuGroups(), src/inner.cpp); no exclusion, no boundary warning.
    theo_sd2 <- nlmixr2data::theo_sd
    theo_sd2$logWT <- log(theo_sd2$WT / 70)

    modBounded <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        allo.cl <- c(0, 0.75, 2)
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

    fitBounded <- .getCachedFit(
      name = "mfocei-bounded-coef",
      fitFn = function() .nlmixr(modBounded, theo_sd2, "mfocei", mfoceiControl(print = 0)),
      cacheFile = "fit-mfocei-bounded-coef.rds"
    )

    # no "cannot respect a boundary" carve-out warning anymore
    expect_false(any(grepl("cannot respect a boundar", fitBounded$runInfo)))
    # the bound is still respected (clamped regression), not silently ignored
    expect_true(unname(fitBounded$theta["allo.cl"]) >= 0)
    expect_true(unname(fitBounded$theta["allo.cl"]) <= 2)
    # regression-updated (profiled out of the outer set)
    expect_true("allo.cl" %in%
                  nlmixr2est:::.foceiMuSkipThetaNames(
                    fitBounded$ui,
                    fitBounded$ui$iniDf$name[!is.na(fitBounded$ui$iniDf$ntheta)]))
  })

  test_that("mfocei mu-references a whole group when a sibling covariate is bounded", {
    theo_sd2 <- nlmixr2data::theo_sd
    theo_sd2$logWT <- log(theo_sd2$WT / 70)
    theo_sd2$sexf <- as.numeric(theo_sd2$ID) %% 2

    modMixed <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        allo.cl <- 0.75 # unbounded
        allo.cl2 <- c(-1, 0.1, 1) # bounded
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + allo.cl * logWT + allo.cl2 * sexf)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    fitMixed <- .getCachedFit(
      name = "mfocei-mixed-bounded-coef",
      fitFn = function() .nlmixr(modMixed, theo_sd2, "mfocei", mfoceiControl(print = 0)),
      cacheFile = "fit-mfocei-mixed-bounded-coef.rds"
    )

    expect_false(any(grepl("cannot respect a boundar", fitMixed$runInfo)))
    expect_true(unname(fitMixed$theta["allo.cl2"]) >= -1)
    expect_true(unname(fitMixed$theta["allo.cl2"]) <= 1)
    # every group parameter gets a real (full-model recompute) standard error
    .pf <- fitMixed$parFixed
    expect_false(is.na(suppressWarnings(as.numeric(.pf["allo.cl", "SE"]))))
    expect_false(is.na(suppressWarnings(as.numeric(.pf["tcl", "SE"]))))
    expect_false(is.na(suppressWarnings(as.numeric(.pf["allo.cl2", "SE"]))))
  })

  test_that("mfocei's live iteration print shows mu-group thetas as standard columns", {
    # The mu-group population/covariate thetas are regression-updated (no
    # outer-optimizer slot) but print as normal scale.h columns at their
    # natural theta positions; gradient rows show blank cells for them.
    #
    # capture.output() must wrap the *un-suppressed* nlmixr() call directly
    # (not the .nlmixr()/suppressMessages() test helper) -- this console
    # output is otherwise swallowed by suppressMessages() in this R session.
    # width=200 keeps all 8 columns on one print row (no wrap).
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

    out <- withr::with_options(list(width = 200), capture.output({
      nlmixr2est::nlmixr(mod, theo_sd2, "mfocei",
                          mfoceiControl(print = 1, maxOuterIterations = 2))
    }))

    # the old bolt-on "|   mu|" row is gone
    expect_false(any(grepl("^\\|   mu\\|", out)))
    # Key legend notes the regression-updated thetas
    expect_true(any(grepl("mu-referenced thetas are regression-updated", out)))
    # mu-group thetas (pop + covariate coefficient) are header columns at
    # their natural positions among the other parameters
    hdr <- grep("^\\|    #\\|", out, value = TRUE)[1]
    expect_false(is.na(hdr))
    .pos <- vapply(c("tka", "tcl", "tv", "allo\\.cl", "add\\.sd"),
                   function(nm) as.numeric(regexpr(paste0("\\b", nm, "\\b"), hdr)),
                   numeric(1))
    expect_true(all(.pos > 0))
    expect_true(all(diff(.pos) > 0))
    # gradient rows: one blank cell per regression-updated theta
    # (tka/tcl/tv/allo.cl), real numbers for the optimizer-owned columns
    gradRows <- grep("^\\|    [GFCMSA]\\|", out, value = TRUE)
    expect_true(length(gradRows) > 0)
    .blanks <- vapply(gradRows,
                      function(r) sum(gregexpr(" {11}\\|", r)[[1]] > 0),
                      numeric(1))
    expect_true(all(.blanks == 4))
    expect_true(any(grepl("[-0-9]\\.[0-9]", gradRows)))
  })
})

nmTest({
  test_that("mfocei handles >=2 mu-ref covariate expressions (#711)", {
    # two mu-referenced covariates that are expressions (log(WT/70)), not
    # bare data columns, used to error with "undefined columns selected"
    mod <- function() {
      ini({
        tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
        wtcl <- 0.75; wtv <- 1
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + wtcl * log(WT/70))
        v <- exp(tv + eta.v + wtv * log(WT/70))
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl/v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .f <- .nlmixr(mod, nlmixr2data::theo_sd, "focei",
                  foceiControl(muModel = "lin", print = 0,
                               maxOuterIterations = 0L, covMethod = "",
                               calcTables = FALSE))
    expect_true(inherits(.f, "nlmixr2FitCore"))
    expect_true(all(is.finite(.f$theta)))
  })
})
