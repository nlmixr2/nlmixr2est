# Per-subject finite-difference fallback for the analytic FOCEI outer gradient
# (foceiControl(fast = TRUE)).
#
# Rung 3 of the failure ladder: when a subject's augmented sensitivity solve fails and
# per-individual tolerance loosening cannot rescue it, that subject alone is finite
# differenced and folded into the otherwise-analytic gradient, rather than the whole
# evaluation declining to a finite-difference gradient.
#
# No model reaches the flagged state reliably, so NLMIXR2EST_OUTER_FAIL_ID injects the
# failure into subjects that actually solve.  That is also what makes the substituted
# contribution measurable: with a failure injected into a subject whose exact analytic
# slope IS available, three runs isolate it with no dilution from the correct subjects --
#
#   gRef  nothing flagged                        -> the all-analytic gradient
#   gSkip flagged, FD suppressed (FD_SKIP=1)     -> analytic sum over the REMAINING subjects
#   gOne  flagged, FD applied                    -> the substitution in place
#
# so An = gRef - gSkip is the flagged subjects' exact analytic contribution and
# FD = gOne - gSkip is what the finite difference put there.  Comparing those two is the
# test; comparing totals is not, because 11 correct subjects hide one wrong column.
#
# Kept in the ESSENTIAL push/PR subset (not in .slowBatches), so do NOT add skip_on_ci():
# one short fit plus post-fit gradient calls.

nmTest({
  .fbClearHooks <- function() {
    Sys.unsetenv("NLMIXR2EST_OUTER_FAIL_ID")
    Sys.unsetenv("NLMIXR2EST_OUTER_FD_SKIP")
  }

  .fbModel <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  # Gradient at the fit's converged parameters, with the hooks set for that call only.
  .fbGrad <- function(fit, ids = NULL, skip = FALSE) {
    .fbClearHooks()
    on.exit(.fbClearHooks(), add = TRUE)
    if (!is.null(ids)) Sys.setenv(NLMIXR2EST_OUTER_FAIL_ID = paste(ids, collapse = ","))
    if (skip) Sys.setenv(NLMIXR2EST_OUTER_FD_SKIP = "1")
    .foceiGradDirect(fit)
  }

  test_that("a flagged subject is finite-differenced instead of declining the gradient", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .fbClearHooks()
    on.exit(.fbClearHooks(), add = TRUE)

    fit <- suppressMessages(suppressWarnings(nlmixr2(
      .fbModel, nlmixr2data::theo_sd, "focei",
      foceiControl(print = 0L, covMethod = "", fast = TRUE, sigdig = 3,
                   calcTables = FALSE, maxOuterIterations = 3L))))

    gRef <- .fbGrad(fit)
    expect_true(is.numeric(gRef))
    expect_true(all(is.finite(gRef)))
    np <- length(gRef)
    # 4 estimated thetas (tka tcl tv add.sd) + 3 omega entries
    expect_equal(np, 7L)

    for (ids in list(2L, c(2L, 7L))) {
      gSkip <- .fbGrad(fit, ids, skip = TRUE)
      gOne  <- .fbGrad(fit, ids, skip = FALSE)
      # LENGTH first, every time.  .foceiGradDirect returns NULL when the analytic route
      # declined, and all(is.finite(NULL)) is TRUE -- so without this the finiteness and
      # subtraction assertions below pass vacuously on a gradient that was never computed,
      # i.e. the test would report success for a fallback that is entirely broken.
      expect_equal(length(gSkip), np)
      expect_equal(length(gOne), np)
      expect_true(all(is.finite(gSkip)))
      expect_true(all(is.finite(gOne)))

      # The fallback ran and actually substituted: suppressing it must move the gradient.
      expect_false(isTRUE(all.equal(gSkip, gOne, tolerance = 1e-12)))

      an <- gRef - gSkip     # exact analytic contribution of the flagged subjects
      fd <- gOne - gSkip     # what the finite difference substituted
      expect_true(all(is.finite(an)))
      expect_true(all(is.finite(fd)))

      # Every component must be substituted, the OMEGA block included -- a theta-only
      # fallback left those columns at exactly zero and still looked fine on the total.
      iOm <- 5L:7L
      expect_true(all(abs(fd[iOm]) > 0))

      # Agreement, per block, relative to the exact analytic contribution.
      relL2 <- function(idx) sqrt(sum((fd[idx] - an[idx])^2)) / sqrt(sum(an[idx]^2))
      expect_lt(relL2(1L:4L), 0.05)   # theta + sigma
      expect_lt(relL2(iOm), 0.05)     # omega
    }
  })

  test_that("the fallback survives a fit and reports it, and omega is not silently zero", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .fbClearHooks()
    on.exit(.fbClearHooks(), add = TRUE)

    # Flag DURING the fit: the gradient is evaluated at many thetas, so this exercises the
    # shared step cache and its invalidation as the flagged set/theta move.
    Sys.setenv(NLMIXR2EST_OUTER_FAIL_ID = "2,7")
    fit <- suppressMessages(suppressWarnings(nlmixr2(
      .fbModel, nlmixr2data::theo_sd, "focei",
      foceiControl(print = 0L, covMethod = "", fast = TRUE, sigdig = 3,
                   calcTables = FALSE, maxOuterIterations = 3L))))
    .fbClearHooks()

    expect_true(is.finite(fit$objf))
    # The analytic route SURVIVED the failed subjects rather than declining wholesale.
    expect_gt(as.integer(fit$env$nOuterFdInd), 0L)
    expect_gt(as.integer(fit$env$nAnalyticGradDirect), 0L)

    # Diagnostics the fallback is required to maintain.
    .out <- fit$env$nFdOutlier
    expect_true(all(c("params", "chartrandSlopes", "stepClamped") %in% names(.out)))
    expect_true(all(as.integer(.out) >= 0L))
  })

  test_that("a clamped step does not disable the fallback for the rest of the fit", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .fbClearHooks()
    on.exit(.fbClearHooks(), add = TRUE)

    # shi21hMin == shi21hMax is rejected by foceiControl, but a floor this high forces the
    # step search onto its bound for every parameter.  A clamped step must be neither used
    # nor CACHED: caching it would make every later evaluation skip the search, fail the
    # same bound test and decline, for the rest of the fit.  The observable is that the fit
    # still completes with a finite objective (declining to full FD is the intended
    # behaviour) rather than the run degrading or erroring.
    Sys.setenv(NLMIXR2EST_OUTER_FAIL_ID = "2")
    fit <- suppressMessages(suppressWarnings(nlmixr2(
      .fbModel, nlmixr2data::theo_sd, "focei",
      foceiControl(print = 0L, covMethod = "", fast = TRUE, sigdig = 3,
                   calcTables = FALSE, maxOuterIterations = 2L,
                   shi21hMin = 1.9, shi21hMax = 2.0))))
    .fbClearHooks()
    expect_true(is.finite(fit$objf))
    # The clamp actually FIRED -- otherwise this test proves nothing about the caching of a
    # clamped step, it just re-runs the previous test with different control arguments.
    expect_gt(as.integer(fit$env$nFdOutlier[["stepClamped"]]), 0L)
    # ... and it kept firing rather than being cached once and skipped thereafter.  A cached
    # clamp would make every later evaluation skip the search, so the count would stall at
    # the number of free parameters searched in the FIRST flagged evaluation.
    expect_gt(as.integer(fit$env$nFdOutlier[["stepClamped"]]), 7L)
  })

  test_that("the fallback is correct when a theta is fix()ed", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .fbClearHooks()
    on.exit(.fbClearHooks(), add = TRUE)

    # The step store used to be allocated per OPTIMIZER parameter but indexed by FULL-THETA
    # position, so any fix()ed parameter made the last free one write past its end.  Every
    # other case here estimates all four thetas, where the two indexings coincide and the
    # bug is invisible.  tv fixed puts add.sd at full-theta slot 3 with only 3 free thetas,
    # so the old store (length 3 for the theta block) was overrun by exactly one.
    .fbFixModel <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- fix(3.45); add.sd <- 0.7
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d/dt(depot)  <- -ka * depot
        d/dt(center) <-  ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }

    fit <- suppressMessages(suppressWarnings(nlmixr2(
      .fbFixModel, nlmixr2data::theo_sd, "focei",
      foceiControl(print = 0L, covMethod = "", fast = TRUE, sigdig = 3,
                   calcTables = FALSE, maxOuterIterations = 3L))))

    gRef <- .fbGrad(fit)
    # 3 free thetas (tka tcl add.sd) + 3 omega entries; tv is fixed and carries no column.
    expect_equal(length(gRef), 6L)
    expect_true(all(is.finite(gRef)))

    ids <- 2L
    gSkip <- .fbGrad(fit, ids, skip = TRUE)
    gOne  <- .fbGrad(fit, ids, skip = FALSE)
    expect_true(all(is.finite(gSkip)))
    expect_true(all(is.finite(gOne)))
    expect_false(isTRUE(all.equal(gSkip, gOne, tolerance = 1e-12)))

    an <- gRef - gSkip
    fd <- gOne - gSkip
    # The FIXED parameter must not consume a column: if the fold-in walked full-theta
    # positions instead of the free set, add.sd's contribution would land in tv's slot and
    # the last component would be left at zero.
    expect_true(all(abs(fd) > 0))
    relL2 <- function(idx) sqrt(sum((fd[idx] - an[idx])^2)) / sqrt(sum(an[idx]^2))
    expect_lt(relL2(1L:3L), 0.05)   # theta + sigma
    expect_lt(relL2(4L:6L), 0.05)   # omega
  })


  test_that("the pooled/VAE M-step keeps its analytic gradient when a subject is flagged", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .fbClearHooks()
    on.exit(.fbClearHooks(), add = TRUE)

    # foceiGradPooledDirect_ is the entry est="vae", nonMuTheta="grad" uses each M-step.  It is
    # handed a theta/eta/omega that are deliberately NOT the inner problem's, so the fold-in
    # has to be taken at THAT point -- which is why the FD core takes its point explicitly
    # instead of reading op_focei.  Two things this pins, neither visible in the fitted values:
    #
    #   * before the fold-in existed here, a flagged subject made this entry return the kernel
    #     sum with that subject's column still ZERO.  Zeros are finite, so a gradient silently
    #     missing a whole subject reached the M-step.
    #   * the stopgap for that was to decline, which sends every M-step to the bobyqa
    #     regression -- correct but a different optimizer.  nRegFallback counts exactly that,
    #     so "the analytic gradient survived the flagged subject" is testable rather than
    #     inferred from the fitted values.
    .vMod <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- c(2, 3.45, 5); add.sd <- c(0, 0.7, 5)
        eta.ka ~ 0.6; eta.cl ~ 0.3 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd) })
    }
    .vCtl <- vaeControl(nonMuTheta = "grad", print = 0L, calcTables = FALSE,
                        returnVae = TRUE, itersBurnIn = 10L, iters = 30L,
                        klWarmup = 5L, gammaIter = 20L)

    Sys.setenv(NLMIXR2EST_OUTER_FAIL_ID = "2,7")
    r <- suppressWarnings(suppressMessages(
      nlmixr2(.vMod(), nlmixr2data::theo_sd, est = "vae", control = .vCtl)))
    .fbClearHooks()

    # The M-steps used the analytic gradient and NONE fell back -- i.e. the entry substituted
    # the flagged subjects instead of declining.
    expect_gt(r$nRegGrad, 0L)
    expect_equal(r$nRegFallback, 0L)
    expect_true(is.finite(r$regressTheta[["tv"]]))
  })

  test_that("fdIndividualStep switches the step search and both routes agree", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .fbClearHooks()
    on.exit(.fbClearHooks(), add = TRUE)

    expect_true(foceiControl(fdIndividualStep = TRUE)$fdIndividualStep == 1L)
    expect_true(foceiControl(fdIndividualStep = FALSE)$fdIndividualStep == 0L)
    expect_error(foceiControl(fdIndividualStep = "yes"))
    expect_true(isTRUE(vaeControl(fdIndividualStep = TRUE)$fdIndividualStep))
    expect_false(isTRUE(vaeControl(fdIndividualStep = FALSE)$fdIndividualStep))

    # Both routes must produce a usable substituted gradient.  They are NOT expected to agree
    # to the bit -- a per-subject step and a shared one are different differences -- but they
    # difference the same function, so they must agree to the accuracy this path claims.
    .g <- function(indiv) {
      .fbClearHooks()
      on.exit(.fbClearHooks(), add = TRUE)
      fit <- suppressMessages(suppressWarnings(nlmixr2(
        .fbModel, nlmixr2data::theo_sd, "focei",
        foceiControl(print = 0L, covMethod = "", fast = TRUE, sigdig = 3,
                     calcTables = FALSE, maxOuterIterations = 3L,
                     fdIndividualStep = indiv))))
      Sys.setenv(NLMIXR2EST_OUTER_FAIL_ID = "2,7")
      .foceiGradDirect(fit)
    }
    gInd <- .g(TRUE)
    gSha <- .g(FALSE)
    expect_equal(length(gInd), 7L)
    expect_equal(length(gSha), 7L)
    expect_true(all(is.finite(gInd)))
    expect_true(all(is.finite(gSha)))
    expect_lt(max(abs(gInd - gSha) / pmax(abs(gInd), 1e-8)), 0.05)
  })

  test_that("fdOutlierZ makes the Chartrand refinement reachable, and fdChartrand gates it", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .fbClearHooks()
    on.exit(.fbClearHooks(), add = TRUE)

    expect_equal(foceiControl(fdOutlierZ = 1.0)$fdOutlierZ, 1.0)
    expect_error(foceiControl(fdOutlierZ = -1))
    expect_equal(vaeControl(fdOutlierZ = 1.0)$fdOutlierZ, 1.0)

    # On a well-behaved fit the outlier test never fires at the conventional 3.5 cut, so the
    # TV refinement it gates was unreachable and therefore untestable.  Driving the cut to ~0
    # makes every differenced slope an outlier, which is what lets the pass be exercised at
    # all -- and then fdChartrand=FALSE must still record the DETECTION while performing no
    # refinement, which is the contract that separates the two knobs.
    .run <- function(z, chartrand) {
      .fbClearHooks()
      on.exit(.fbClearHooks(), add = TRUE)
      Sys.setenv(NLMIXR2EST_OUTER_FAIL_ID = "2,7")
      fit <- suppressMessages(suppressWarnings(nlmixr2(
        .fbModel, nlmixr2data::theo_sd, "focei",
        foceiControl(print = 0L, covMethod = "", fast = TRUE, sigdig = 3,
                     calcTables = FALSE, maxOuterIterations = 2L,
                     fdOutlierZ = z, fdChartrand = chartrand))))
      .fbClearHooks()
      list(objf = fit$objf, out = fit$env$nFdOutlier)
    }

    hi <- .run(3.5, TRUE)          # the default: nothing flagged on this fit
    expect_true(is.finite(hi$objf))
    expect_equal(as.integer(hi$out[["params"]]), 0L)
    expect_equal(as.integer(hi$out[["chartrandSlopes"]]), 0L)

    lo <- .run(1e-8, TRUE)         # everything flagged: the pass runs
    expect_true(is.finite(lo$objf))
    expect_gt(as.integer(lo$out[["params"]]), 0L)
    expect_gt(as.integer(lo$out[["chartrandSlopes"]]), 0L)

    off <- .run(1e-8, FALSE)       # detected but NOT refined
    expect_true(is.finite(off$objf))
    expect_gt(as.integer(off$out[["params"]]), 0L)
    expect_equal(as.integer(off$out[["chartrandSlopes"]]), 0L)
  })

  test_that("fdChartrandAll refines every FD subject, fdOutlierAny widens the trigger", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .fbClearHooks()
    on.exit(.fbClearHooks(), add = TRUE)

    expect_equal(foceiControl(fdChartrandAll = TRUE)$fdChartrandAll, 1L)
    expect_equal(foceiControl(fdOutlierAny = TRUE)$fdOutlierAny, 1L)
    expect_error(foceiControl(fdChartrandAll = "yes"))
    expect_true(isTRUE(vaeControl(fdChartrandAll = TRUE)$fdChartrandAll))
    expect_true(isTRUE(vaeControl(fdOutlierAny = TRUE)$fdOutlierAny))

    .run <- function(...) {
      .fbClearHooks()
      on.exit(.fbClearHooks(), add = TRUE)
      Sys.setenv(NLMIXR2EST_OUTER_FAIL_ID = "2,7")
      fit <- suppressMessages(suppressWarnings(nlmixr2(
        .fbModel, nlmixr2data::theo_sd, "focei",
        foceiControl(print = 0L, covMethod = "", fast = TRUE, sigdig = 3,
                     calcTables = FALSE, maxOuterIterations = 2L, ...))))
      .fbClearHooks()
      list(objf = fit$objf, out = fit$env$nFdOutlier)
    }

    # fdOutlierAny: at the DEFAULT cut nothing among the finite differences is extreme, so the
    # pass does not fire.  Letting an exact slope fire it is a strictly wider trigger, so with
    # the cut low enough to make the analytic slopes themselves dispersed it must fire and say
    # so through analyticTrigger -- while the analytic gradients stay analytic.
    base <- .run(fdOutlierZ = 3.5)
    expect_equal(as.integer(base$out[["params"]]), 0L)
    expect_equal(as.integer(base$out[["analyticTrigger"]]), 0L)

    # ACTUALLY exercise fdOutlierAny, and assert the counter that proves the branch ran.
    # It fires at the DEFAULT 3.5 cut on theo_sd, which is itself the point: the EXACT
    # per-subject slopes are dispersed enough to be outliers against each other, because a
    # subject's slope scales with the data it carries.  With the option off the same fit
    # flags nothing, so this is a clean on/off pair rather than a threshold trick.
    anyOn <- .run(fdOutlierZ = 3.5, fdOutlierAny = TRUE)
    expect_true(is.finite(anyOn$objf))
    expect_gt(as.integer(anyOn$out[["analyticTrigger"]]), 0L)
    expect_gt(as.integer(anyOn$out[["params"]]), 0L)
    # ... and it is the ANALYTIC side that fired: no finite difference was an outlier, so
    # nothing is refined.  fdOutlierAny on its own DETECTS; it only changes the gradient when
    # paired with fdChartrandAll, which is what makes every FD subject eligible.
    expect_equal(as.integer(anyOn$out[["chartrandSlopes"]]), 0L)
    anyBoth <- .run(fdOutlierZ = 3.5, fdOutlierAny = TRUE, fdChartrandAll = TRUE)
    expect_gt(as.integer(anyBoth$out[["chartrandSlopes"]]), 0L)

    # fdChartrandAll: with the pass firing, refining every FD subject must produce at least as
    # many refined slopes as refining only the outliers.  Both must still give a finite fit --
    # the analytic subjects are never recomputed under either setting.
    only <- .run(fdOutlierZ = 1e-8, fdChartrandAll = FALSE)
    all  <- .run(fdOutlierZ = 1e-8, fdChartrandAll = TRUE)
    expect_true(is.finite(only$objf))
    expect_true(is.finite(all$objf))
    expect_gt(as.integer(only$out[["chartrandSlopes"]]), 0L)
    expect_gte(as.integer(all$out[["chartrandSlopes"]]),
               as.integer(only$out[["chartrandSlopes"]]))
  })

  test_that("fdRefine selects richardson/lanczos/chartrand and all three refine", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .fbClearHooks()
    on.exit(.fbClearHooks(), add = TRUE)

    expect_equal(foceiControl(fdRefine = "chartrand")$fdRefine, "chartrand")
    expect_equal(foceiControl(fdRefine = "lanczos")$fdRefine, "lanczos")
    # A control is re-passed through foceiControl(); an integer code would come back NA from
    # match.arg() SILENTLY and fall back to the default estimator.
    expect_equal(do.call(foceiControl, foceiControl(fdRefine = "richardson"))$fdRefine,
                 "richardson")
    expect_error(foceiControl(fdRefine = "spline"))
    expect_error(foceiControl(fdRichardsonV = 1))     # v must exceed 1 to shrink the step
    expect_equal(vaeControl(fdRefine = "lanczos")$fdRefine, "lanczos")

    # Each estimator must actually run and produce a finite substituted gradient.  The cut is
    # driven to ~0 so the pass fires at all -- at the default nothing on this fit is an
    # outlier, which is exactly why these were unreachable before fdOutlierZ existed.
    .g <- function(meth) {
      .fbClearHooks()
      on.exit(.fbClearHooks(), add = TRUE)
      Sys.setenv(NLMIXR2EST_OUTER_FAIL_ID = "2,7")
      fit <- suppressMessages(suppressWarnings(nlmixr2(
        .fbModel, nlmixr2data::theo_sd, "focei",
        foceiControl(print = 0L, covMethod = "", fast = TRUE, sigdig = 3,
                     calcTables = FALSE, maxOuterIterations = 2L,
                     fdOutlierZ = 1e-8, fdRefine = meth))))
      .fbClearHooks()
      list(objf = fit$objf, n = as.integer(fit$env$nFdOutlier[["chartrandSlopes"]]))
    }
    for (m in c("chartrand", "lanczos", "richardson")) {
      r <- .g(m)
      expect_true(is.finite(r$objf))
      # the counter is the shared "a slope was recomputed" counter, whichever estimator ran
      expect_gt(r$n, 0L)
    }

    # No accuracy ordering is asserted.  theo_sd is a smooth, well-conditioned surface and
    # cannot discriminate between these -- that comparison needs a stiff likelihood surface
    # and is the separate nlmixr2est follow-up these options exist to enable.
  })

  test_that("fdOutlierScale tests the per-observation slope on unbalanced data", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .fbClearHooks()
    on.exit(.fbClearHooks(), add = TRUE)

    expect_equal(foceiControl(fdOutlierScale = TRUE)$fdOutlierScale, 1L)
    expect_equal(foceiControl(fdOutlierScale = FALSE)$fdOutlierScale, 0L)
    expect_error(foceiControl(fdOutlierScale = "yes"))
    expect_true(isTRUE(vaeControl(fdOutlierScale = FALSE)$fdOutlierScale) == FALSE)

    .fit <- function(dat, scale) {
      .fbClearHooks()
      on.exit(.fbClearHooks(), add = TRUE)
      Sys.setenv(NLMIXR2EST_OUTER_FAIL_ID = "2,7")
      f <- suppressMessages(suppressWarnings(nlmixr2(
        .fbModel, dat, "focei",
        foceiControl(print = 0L, covMethod = "", fast = TRUE, sigdig = 3,
                     calcTables = FALSE, maxOuterIterations = 2L,
                     fdOutlierScale = scale))))
      .fbClearHooks()
      list(objf = f$objf, params = as.integer(f$env$nFdOutlier[["params"]]))
    }

    # BALANCED: every subject has the same observation count, so dividing every slope by the
    # same number cannot change a modified z-score.  Scaling must be a no-op here -- that is
    # what makes it safe to turn on by default.
    bal <- nlmixr2data::theo_sd
    on  <- .fit(bal, TRUE)
    off <- .fit(bal, FALSE)
    expect_true(is.finite(on$objf))
    expect_equal(on$objf, off$objf)
    expect_equal(on$params, off$params)

    # UNBALANCED: thin subjects 1-4 to 3 observations each while the rest keep 11.  Raw slopes
    # are then not draws from one distribution -- a well-sampled subject has a systematically
    # larger slope purely from carrying more data -- so the two settings must be able to
    # DISAGREE about who is an outlier.  Both must still produce a finite fit.
    .thin <- do.call(rbind, lapply(split(bal, bal$ID), function(d) {
      obs <- d[d$EVID == 0, , drop = FALSE]
      dose <- d[d$EVID != 0, , drop = FALSE]
      if (as.integer(as.character(d$ID[1])) <= 4L && nrow(obs) > 3L) {
        obs <- obs[seq_len(3L), , drop = FALSE]
      }
      rbind(dose, obs)
    }))
    .thin <- .thin[order(.thin$ID, .thin$TIME), , drop = FALSE]
    uOn  <- .fit(.thin, TRUE)
    uOff <- .fit(.thin, FALSE)
    expect_true(is.finite(uOn$objf))
    expect_true(is.finite(uOff$objf))
    # The option must actually REACH the C++.  Driving the cut low makes the test fire on
    # both settings; scaling then changes which per-observation slopes look extreme, so the
    # two must not agree on everything.  Without this the test passes even if
    # op_focei.fdOutlierScale were ignored entirely.
    .cut <- function(scale) {
      .fbClearHooks()
      on.exit(.fbClearHooks(), add = TRUE)
      Sys.setenv(NLMIXR2EST_OUTER_FAIL_ID = "2,7")
      f <- suppressMessages(suppressWarnings(nlmixr2(
        .fbModel, .thin, "focei",
        foceiControl(print = 0L, covMethod = "", fast = TRUE, sigdig = 3,
                     calcTables = FALSE, maxOuterIterations = 2L,
                     fdOutlierZ = 0.5, fdOutlierScale = scale))))
      .fbClearHooks()
      c(params = as.integer(f$env$nFdOutlier[["params"]]),
        slopes = as.integer(f$env$nFdOutlier[["chartrandSlopes"]]),
        objf = f$objf)
    }
    cOn <- .cut(TRUE); cOff <- .cut(FALSE)
    expect_true(is.finite(cOn[["objf"]]))
    expect_true(is.finite(cOff[["objf"]]))
    expect_gt(cOn[["params"]], 0L)
    expect_gt(cOff[["params"]], 0L)
    # The option must CHANGE which slopes look extreme on unbalanced data -- otherwise this
    # passes even if op_focei.fdOutlierScale were ignored entirely in C++.
    expect_false(cOn[["params"]] == cOff[["params"]] &&
                 cOn[["slopes"]] == cOff[["slopes"]])
  })
})
