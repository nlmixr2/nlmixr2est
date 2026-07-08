nmTest({
  test_that("mu/irls covariance matches the full base focei model (linCmt/ODE/matExp)", {
    # The mu/irls covariance step bails and is recomputed on the full base model
    # (muModel="none"); the reported SEs -- including the mu-ref/covariate ("linear")
    # parameters -- must match the equivalent plain focei fit, not the mu->phi reduction.
    theo_sd2 <- nlmixr2data::theo_sd
    theo_sd2$logWT <- log(theo_sd2$WT / 70)

    modLin <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; allo.cl <- 0.75
            eta.ka ~ 0.6; eta.cl ~ 0.3; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl + allo.cl * logWT); v <- exp(tv)
        linCmt() ~ add(add.sd) })
    }
    modOde <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; allo.cl <- 0.75
            eta.ka ~ 0.6; eta.cl ~ 0.3; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl + allo.cl * logWT); v <- exp(tv)
        d/dt(depot) <- -ka * depot; d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v; cp ~ add(add.sd) })
    }
    # non-covariate matExp: exercises the full corresponding model on the matrix-
    # exponential path (the covariate coefficient is dropped so the comparison is clean)
    matNoCov <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; eta.cl ~ 0.3; add.sd <- 0.7 })
      model({ matExp()
        k_depot_central <- exp(tka + eta.ka)
        k_central_output <- exp(tcl + eta.cl) / exp(tv)
        cp <- central / exp(tv); cp ~ add(add.sd) })
    }

    .chk <- function(mfun, muEst, muCtl, hasCov = TRUE) {
      .fB <- .nlmixr(mfun, theo_sd2, "focei", foceiControl(print = 0, covMethod = "r,s"))
      .fM <- .nlmixr(mfun, theo_sd2, muEst, muCtl(print = 0, covMethod = "r,s"))
      .sB <- sqrt(diag(.fB$cov)); .sM <- sqrt(diag(.fM$cov))
      # the mu fit takes a real covariance (not the bailed empty one): finite SEs
      expect_false(is.na(suppressWarnings(as.numeric(.fM$parFixed["tcl", "SE"]))))
      # well-identified structural thetas match the plain focei fit closely.  tcl was
      # the ~8.5x-wrong parameter under the mu->phi reduction; it now matches.
      .th <- intersect(c("tka", "tcl", "tv", "add.sd"), names(.sB))
      expect_equal(unname(.sM[.th]), unname(.sB[.th]), tolerance = 0.1)
      if (hasCov) {
        # covariate ("linear") coefficient: finite SE and the right order of magnitude
        # (the bug made it ~8x off; focei's optimizer and the mu regression converge the
        # weakly-identified allo.cl to slightly different points, so allow a factor of 2)
        expect_false(is.na(suppressWarnings(as.numeric(.fM$parFixed["allo.cl", "SE"]))))
        .r <- as.numeric(.sM["allo.cl"]) / as.numeric(.sB["allo.cl"])
        expect_true(is.finite(.r) && .r > 0.5 && .r < 2)
      }
    }

    for (.est in c("mufocei", "irlsfocei")) {
      .ctl <- if (.est == "mufocei") mufoceiControl else irlsfoceiControl
      .chk(modLin, .est, .ctl)                 # linCmt covariate
      .chk(modOde, .est, .ctl)                 # ODE covariate
      .chk(matNoCov, .est, .ctl, hasCov = FALSE)  # matExp -- full corresponding model
    }
  })

  test_that("FOCE/focep mu families recompute covariance on the full base model (ODE/matExp)", {
    # The FOCE- and foce+-based mu/irls variants pivot to foce/focep for the cov step.
    # Use well-identified (no-covariate) models so foce/focep and the mu regression
    # converge to the same point; the covariate case is exercised for FOCEI above.
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

    .chk <- function(mfun, muEst, baseEst, muCtl, baseCtl) {
      .fB <- .nlmixr(mfun, theo_sd2, baseEst, baseCtl(print = 0, covMethod = "r,s"))
      .fM <- .nlmixr(mfun, theo_sd2, muEst, muCtl(print = 0, covMethod = "r,s"))
      expect_false(is.na(suppressWarnings(as.numeric(.fM$parFixed["tcl", "SE"]))))
      .sB <- sqrt(diag(.fB$cov)); .sM <- sqrt(diag(.fM$cov))
      .cmn <- intersect(names(.sB), names(.sM))
      expect_equal(unname(.sM[.cmn]), unname(.sB[.cmn]), tolerance = 0.05)
    }

    for (.mod in list(odeNC, matNC)) {
      .chk(.mod, "mufoce",    "foce",  mufoceControl,    foceControl)
      .chk(.mod, "irlsfoce",  "foce",  irlsfoceControl,  foceControl)
      .chk(.mod, "mufocep",   "focep", mufocepControl,   focepControl)
      .chk(.mod, "irlsfocep", "focep", irlsfocepControl, focepControl)
    }
  })

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

  test_that("mufocei excludes a bounded mu-ref covariate coefficient, warns, and still respects the bound", {
    # The closed-form/IRLS regression cannot respect a box constraint on
    # a covariate coefficient, so it is treated as if it were
    # time-varying and carved out of the regression -- estimated as an
    # ordinary bounded theta by the outer optimizer instead
    # (.muRefGroups(), R/muRefClassify.R) -- with a warning explaining why,
    # captured into the fit's runInfo the same way every other pre-fit
    # warning is (R/nlmixr2Est.R's .collectWarn() wraps the whole
    # nlmixr2Est() dispatch). The group's population theta still benefits
    # from the mu-ref speed-up (only the bounded slope is excluded).
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
      name = "mufocei-bounded-coef",
      fitFn = function() .nlmixr(modBounded, theo_sd2, "mufocei", mufoceiControl(print = 0)),
      cacheFile = "fit-mufocei-bounded-coef.rds"
    )

    expect_true(any(grepl("allo\\.cl.*boundar", fitBounded$runInfo)))
    # the bound is still respected (ordinary bounded-optimizer handling),
    # not silently ignored
    expect_true(unname(fitBounded$theta["allo.cl"]) >= 0)
    expect_true(unname(fitBounded$theta["allo.cl"]) <= 2)
  })

  test_that("mufocei keeps mu-referencing a group's unbounded covariate/population theta when only a sibling covariate is bounded", {
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
      name = "mufocei-mixed-bounded-coef",
      fitFn = function() .nlmixr(modMixed, theo_sd2, "mufocei", mufoceiControl(print = 0)),
      cacheFile = "fit-mufocei-mixed-bounded-coef.rds"
    )

    expect_true(any(grepl("allo\\.cl2.*boundar", fitMixed$runInfo)))
    expect_true(unname(fitMixed$theta["allo.cl2"]) >= -1)
    expect_true(unname(fitMixed$theta["allo.cl2"]) <= 1)
    # the unbounded covariate and the group's population theta still get
    # a real (mu-ref-derived) standard error, unaffected by the sibling
    # covariate's exclusion
    .pf <- fitMixed$parFixed
    expect_false(is.na(suppressWarnings(as.numeric(.pf["allo.cl", "SE"]))))
    expect_false(is.na(suppressWarnings(as.numeric(.pf["tcl", "SE"]))))
    # the bounded covariate is an ordinary theta -- it gets a ordinary SE too
    expect_false(is.na(suppressWarnings(as.numeric(.pf["allo.cl2", "SE"]))))
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
