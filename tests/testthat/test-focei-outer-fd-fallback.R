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
})
