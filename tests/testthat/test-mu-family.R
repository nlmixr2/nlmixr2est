nmTest({
  test_that("mu-referenced FOCEI family: all 7 generalized methods recover comparable estimates to their plain counterparts", {
    # Local model/fits (not the shared helper-zzz-fits.R cache): these tests
    # are specifically about the mu-referenced FOCEI-family mechanism
    # itself, per the "WHEN TO KEEP LOCAL FITS" guidance in
    # helper-zzz-fits.R. Same theo_sd-based covariate-on-clearance model
    # used by test-mfocei.R, reused here so all 8 family members and
    # their 4 plain counterparts are validated on identical model/data.

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

    .fitPlain <- function(est, ctl) {
      .getCachedFit(
        name = paste0("mu-family-", est),
        fitFn = function() .nlmixr(mod, theo_sd2, est, ctl),
        cacheFile = paste0("fit-mu-family-", est, ".rds")
      )
    }

    # base ("plain") fits shared across the comparisons below
    fitFocei <- .fitPlain("focei", foceiControl(print = 0))
    fitFoce <- .fitPlain("foce", foceControl(print = 0))
    fitAgq <- .fitPlain("agq", agqControl(print = 0, nAGQ = 2))
    fitLaplace <- .fitPlain("laplace", laplaceControl(print = 0))

    # one (plainFit, muFit, tolerance-on-allo.cl) triple per new method;
    # the other theta/omega/objf tolerances match test-mfocei.R's
    # already-validated values
    .cases <- list(
      ifocei = list(plain = fitFocei, ctl = ifoceiControl(print = 0)),
      mfoce = list(plain = fitFoce, ctl = mfoceControl(print = 0)),
      ifoce = list(plain = fitFoce, ctl = ifoceControl(print = 0)),
      magq = list(plain = fitAgq, ctl = magqControl(print = 0, nAGQ = 2)),
      iagq = list(plain = fitAgq, ctl = iagqControl(print = 0, nAGQ = 2)),
      mlaplace = list(plain = fitLaplace, ctl = mlaplaceControl(print = 0)),
      ilaplace = list(plain = fitLaplace, ctl = ilaplaceControl(print = 0))
    )

    for (.est in names(.cases)) {
      .c <- .cases[[.est]]
      .fitMu <- .getCachedFit(
        name = paste0("mu-family-", .est),
        fitFn = function() .nlmixr(mod, theo_sd2, .est, .c$ctl),
        cacheFile = paste0("fit-mu-family-", .est, ".rds")
      )

      # structural: a normal, fully-classed FOCEI-family fit object
      expect_true(inherits(.fitMu, "nlmixr2FitData"), info = .est)
      expect_true(inherits(.fitMu, "nlmixr2FitCore"), info = .est)

      # the mu-ref-covariate theta's coefficient gets a real (non-NA)
      # standard error from the Hessian-unlock step (Phase 4)
      .pf <- .fitMu$parFixed
      expect_true("allo.cl" %in% rownames(.pf), info = .est)
      expect_false(is.na(suppressWarnings(as.numeric(.pf["allo.cl", "SE"]))), info = .est)

      .thPlain <- .c$plain$theta
      .thMu <- .fitMu$theta
      expect_equal(unname(.thMu["tka"]), unname(.thPlain["tka"]), tolerance = 0.05, info = .est)
      expect_equal(unname(.thMu["tcl"]), unname(.thPlain["tcl"]), tolerance = 0.05, info = .est)
      expect_equal(unname(.thMu["tv"]), unname(.thPlain["tv"]), tolerance = 0.05, info = .est)
      expect_equal(unname(.thMu["allo.cl"]), unname(.thPlain["allo.cl"]), tolerance = 0.25, info = .est)
      expect_equal(unname(.thMu["add.sd"]), unname(.thPlain["add.sd"]), tolerance = 0.05, info = .est)

      .omPlain <- diag(.c$plain$omega)
      .omMu <- diag(.fitMu$omega)
      expect_equal(.omMu, .omPlain, tolerance = 0.25, info = .est)

      # objective function values should be in the same ballpark (not an
      # exact match -- different optimization paths)
      expect_equal(.fitMu$objf, .c$plain$objf, tolerance = 0.1, info = .est)
    }
  })

  test_that("mu-referenced FOCEI family: muModel='none' default path is unaffected", {
    # Guards against the Phase 2-4 C++ changes (muGroupThetaSkip,
    # updateMuGroups(), the foceiCalcCov() Hessian-unlock) leaking into the
    # default (muModel="none") path used by every pre-existing FOCEI-family
    # method. Two checks:

    # (a) a model with NO mu-ref covariates at all (one.compartment, the
    # most widely shared fixture in the test suite) still reproduces the
    # long-lived cached fixture computed independently of this feature --
    # if a Phase 2-4 change had altered default behavior even slightly,
    # this would drift.
    fitFresh <- .nlmixr(one.compartment, theo_sd, est = "focei",
                        control = foceiControl(print = 0, maxOuterIterations = 0L))
    expect_equal(unname(fitFresh$theta), unname(one.compartment.fit.focei$theta),
                 tolerance = 1e-6)
    expect_equal(fitFresh$objf, one.compartment.fit.focei$objf, tolerance = 1e-6)

    # (b) a model that DOES have a mu-ref covariate relationship, fit with
    # est="focei" and an EXPLICIT muModel="none", must exactly match the
    # same model/data fit with the implicit default (muModel="none" is
    # foceiControl()'s own default) -- proves the mu-ref-covariate
    # *detection* machinery (.muRefCppGroupSetup(), rxUiGet.foceiOptEnv's
    # foceiMuModel/foceiMuGroup* wiring) is a complete no-op when the
    # feature isn't requested, even on a model where it *could* apply.
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
    fitDefault <- .getCachedFit(
      name = "mu-family-none-default",
      fitFn = function() {
        .nlmixr(mod, theo_sd2, "focei",
                foceiControl(print = 0, maxOuterIterations = 0L))
      },
      cacheFile = "fit-mu-family-none-default.rds"
    )
    fitExplicitNone <- .getCachedFit(
      name = "mu-family-none-explicit",
      fitFn = function() {
        .nlmixr(mod, theo_sd2, "focei",
                foceiControl(print = 0, maxOuterIterations = 0L, muModel = "none"))
      },
      cacheFile = "fit-mu-family-none-explicit.rds"
    )
    expect_equal(unname(fitExplicitNone$theta), unname(fitDefault$theta), tolerance = 1e-8)
    expect_equal(fitExplicitNone$objf, fitDefault$objf, tolerance = 1e-8)
  })

  test_that("mu-referenced FOCEI family: irls remains sane under unequal per-subject information", {
    # Phase 0's original validation goal for "irls" was accuracy at least
    # as good as "lin" under heteroscedastic per-subject information. A
    # strict superiority assertion (irls error < lin error) is too fragile
    # for a regression suite -- nonlinear-optimizer outcomes on a single
    # simulated dataset can go either way by chance even when the
    # underlying method is sound. A ground-truth-recovery assertion is
    # also too fragile here -- with only 24 subjects, half contributing a
    # single observation, the simulated dataset itself is not necessarily
    # well-identified, so both methods can legitimately land on the same
    # small-sample-biased estimate. Instead this checks against a plain
    # `focei` fit of the *same* simulated data as the reliable reference
    # (already extensively validated above and in test-mfocei.R to find
    # the joint-likelihood optimum): "lin" and "irls" should both land
    # close to focei's own answer, i.e. irls's reweighting doesn't
    # destabilize the estimate relative to plain OLS when the per-subject
    # information is uneven.
    .testSeed(1101)
    nsub <- 24
    ids <- seq_len(nsub)
    logWT <- rnorm(nsub, 0, 0.3)
    # half the subjects get 1 observation (low information for their eta),
    # half get 8 (high information)
    nobs <- ifelse(ids <= nsub / 2, 1L, 8L)
    d <- do.call(rbind, lapply(ids, function(i) {
      times <- sort(c(0.25, seq_len(nobs[i] - 1) * 1.5 + 0.25))[seq_len(nobs[i])]
      rbind(
        data.frame(ID = i, TIME = 0, AMT = 320, DV = 0, EVID = 101, logWT = logWT[i]),
        data.frame(ID = i, TIME = times, AMT = 0, DV = 0, EVID = 0, logWT = logWT[i])
      )
    }))

    simMod <- function() {
      ini({
        tka <- 0.5
        tcl <- 1
        tv <- 3.4
        allo.cl <- 0.75
        eta.ka ~ 0.4
        eta.cl ~ 0.2
        eta.v ~ 0.1
        add.sd <- 0.3
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + allo.cl * logWT)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    simUi <- rxode2::rxode2(simMod)
    # rxSolve() on an nlmixr ui with an add(add.sd) endpoint samples both
    # the between-subject etas and the residual error in one pass (`sim`);
    # do not add extra noise on top of it.
    simDat <- rxode2::rxSolve(simUi, d, addDosing = TRUE)
    d$DV <- simDat$sim
    d$DV[d$EVID != 0] <- 0

    fitRef <- .getCachedFit(
      name = "mu-family-hetero-focei",
      fitFn = function() .nlmixr(simMod, d, "focei", foceiControl(print = 0)),
      cacheFile = "fit-mu-family-hetero-focei.rds"
    )
    fitLin <- .getCachedFit(
      name = "mu-family-hetero-lin",
      fitFn = function() .nlmixr(simMod, d, "mfocei", mfoceiControl(print = 0)),
      cacheFile = "fit-mu-family-hetero-lin.rds"
    )
    fitIrls <- .getCachedFit(
      name = "mu-family-hetero-irls",
      fitFn = function() .nlmixr(simMod, d, "ifocei", ifoceiControl(print = 0)),
      cacheFile = "fit-mu-family-hetero-irls.rds"
    )

    .refAllo <- unname(fitRef$theta["allo.cl"])
    expect_equal(unname(fitLin$theta["allo.cl"]), .refAllo, tolerance = 0.3)
    expect_equal(unname(fitIrls$theta["allo.cl"]), .refAllo, tolerance = 0.3)
  })
})
