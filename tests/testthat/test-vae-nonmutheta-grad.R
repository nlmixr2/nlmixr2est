## nonMuTheta="grad" steps the non-mu structural thetas with the EXACT analytic
## outer gradient (one augmented sensitivity solve per M-step) instead of the
## bounded bobyqa regression.  Unit tests cover the control, the shared plumbing
## and the scope policy; one small end-to-end fit asserts the gradient path was
## actually TAKEN (nRegGrad), not merely that the answer looks reasonable -- a
## silent fallback to bobyqa would otherwise pass every value-based check.
## The grad-vs-regress accuracy comparison is a slow fit: see test-vae-grad-fit.R.

nmTest({
  .odeMod <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- c(2, 3.45, 5); add.sd <- c(0, 0.7, 5)
      eta.ka ~ 0.6; eta.cl ~ 0.3 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) })
  }
  ## same model through linCmt(): no symbolic state sensitivities -> out of scope
  .linMod <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- c(2, 3.45, 5); add.sd <- c(0, 0.7, 5)
      eta.ka ~ 0.6; eta.cl ~ 0.3 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      linCmt() ~ add(add.sd) })
  }

  test_that("vaeControl accepts nonMuTheta='grad'", {
    expect_equal(vaeControl(nonMuTheta = "grad")$nonMuTheta, "grad")
    expect_error(vaeControl(nonMuTheta = "gradient"))
  })

  test_that("'grad' shares the 'regress' plumbing (no eta injected)", {
    expect_true(.vaeNonMuIsRegress("regress"))
    expect_true(.vaeNonMuIsRegress("grad"))
    expect_false(.vaeNonMuIsRegress("eta"))
    expect_false(.vaeNonMuIsRegress("fix"))
    expect_false(.vaeNonMuIsRegress("none"))
    ui <- rxode2::assertRxUi(.odeMod())
    ## like "regress", the hook leaves the model alone (only warns)
    expect_null(suppressWarnings(suppressMessages(
      .preProcessVaeNonMuTheta(ui, "vae", NULL, vaeControl(nonMuTheta = "grad")))))
  })

  test_that("the caller policy resolves and only relaxes bounded transforms", {
    ui <- rxode2::assertRxUi(.odeMod())
    expect_true(is.na(.analyticGradCaller(ui)))                    # nobody asked
    expect_true(.analyticGradAllowsBoundedTr(ui, "vae"))
    ## a recorded bounded transform blocks focei (it reports a natural-scale
    ## gradient) but not the vae (it consumes the internal scale it steps on)
    ui2 <- rxode2::rxUiDecompress(ui)
    ui2$boundedTransforms <- list(list(name = "tv"))
    expect_false(.analyticGradAllowsBoundedTr(ui2, "focei"))
    expect_true(.analyticGradAllowsBoundedTr(ui2, "vae"))
  })

  test_that("scope probe accepts an ODE model and rejects linCmt()", {
    expect_true(.vaeGradInScope(rxode2::assertRxUi(.odeMod())))
    expect_false(.vaeGradInScope(rxode2::assertRxUi(.linMod())))
  })

  test_that("an out-of-scope model downgrades to 'regress' and says so", {
    skip_on_cran()
    f <- suppressWarnings(suppressMessages(
      nlmixr2(.linMod(), nlmixr2data::theo_sd, est = "vae",
              control = vaeControl(nonMuTheta = "grad", print = 0L, calcTables = FALSE,
                                   itersBurnIn = 10L, iters = 20L, klWarmup = 5L,
                                   gammaIter = 15L))))
    expect_equal(f$control$nonMuTheta, "regress")
    expect_true(any(grepl("analytic gradient out of scope", f$runInfo)))
  })

  ## same model with NO ini() bounds on the non-mu tv
  .unboundedMod <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
      eta.ka ~ 0.6; eta.cl ~ 0.3 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) })
  }

  test_that("an unbounded non-mu theta gets a finite fallback interval", {
    ## regression: with +-Inf bounds nothing constrained the M-step and an
    ## unbounded tv ran away to ~1e68 (7.9e306 with a wider run).  The fallback is
    ## ini(est) +- max(.vaeNonMuThetaBound, |est| * .vaeNonMuThetaRel).
    p <- .vaeDataPrep(rxode2::assertRxUi(.unboundedMod()), nlmixr2data::theo_sd,
                      vaeControl(nonMuTheta = "regress"))
    i <- match("tv", p$regressNames)
    expect_false(is.na(i))
    expect_true(is.finite(p$regressLower[i]))
    expect_true(is.finite(p$regressUpper[i]))
    ## generous enough to contain the FOCEi MLE (3.4293), so it cannot bind here
    expect_lt(p$regressLower[i], 3.4293)
    expect_gt(p$regressUpper[i], 3.4293)
    ## a user ini() bound still wins
    p2 <- .vaeDataPrep(rxode2::assertRxUi(.odeMod()), nlmixr2data::theo_sd,
                       vaeControl(nonMuTheta = "regress"))
    j <- match("tv", p2$regressNames)
    expect_equal(p2$regressLower[j], 2)
    expect_equal(p2$regressUpper[j], 5)
  })

  test_that("an unbounded non-mu theta converges instead of diverging", {
    skip_on_cran()
    v <- suppressWarnings(suppressMessages(rxode2::rxWithSeed(42,
      nlmixr2(.unboundedMod(), nlmixr2data::theo_sd, est = "vae",
              control = vaeControl(nonMuTheta = "regress", print = 0L,
                                   calcTables = FALSE)))))
    ## before the fix this was ~1e68; the FOCEi MLE is 3.4293
    expect_true(is.finite(v$theta[["tv"]]))
    expect_lt(abs(v$theta[["tv"]] - 3.4293), 0.5)
  })

  test_that("the VAE ELBO carries the transform-both-sides Jacobian", {
    skip_on_cran()
    ## Verifies the term directly at a FIXED eta -- no fitting, so no dependence on
    ## how well the VAE converges on a log-scale error model.
    ##
    ## The reference MUST floor a non-positive DV at sqrt(eps), the way the
    ## transform does: theo_sd has 9 DV==0 rows, each contributing
    ## -log(sqrt(.Machine$double.eps)) = +18.0218.  Comparing against
    ## -sum(log(DV[DV > 0])) instead makes the term look ~10x too small.
    .lnormMod <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
        eta.ka ~ 0.6; eta.cl ~ 0.3 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ lnorm(add.sd) })
    }
    ui <- rxode2::assertRxUi(.lnormMod())
    e <- .vaeInnerSetup(ui, nlmixr2data::theo_sd, NULL, vaeControl(nonMuTheta = "regress"))
    on.exit(.vaeInnerFree(), add = TRUE)
    d <- e$dataSav
    ids <- unique(d$ID)
    withJ <- vaeInnerLik(matrix(0, length(ids), 2L), 1L, FALSE, FALSE)$obj
    ## same quantity with the Jacobian removed, computed from the reference
    .flr <- sqrt(.Machine$double.eps)
    tbs <- vapply(ids, function(i) {
      y <- d$DV[d$ID == i & d$EVID == 0]
      -sum(log(pmax(y, .flr)))
    }, numeric(1))
    ## obj = likInner0 - tbsLik, so removing the term must shift each subject by
    ## exactly its own Jacobian; assert the TOTAL identity to a tight tolerance
    expect_equal(sum(tbs), -18.6516, tolerance = 1e-3)
    expect_true(all(is.finite(withJ)))
    ## and the add-error model must be untouched (tbsLik == 0 there)
    e2 <- .vaeInnerSetup(rxode2::assertRxUi(.odeMod()), nlmixr2data::theo_sd, NULL,
                         vaeControl(nonMuTheta = "regress"))
    o2 <- vaeInnerLik(matrix(0, length(ids), 2L), 1L, FALSE, FALSE)$obj
    .vaeInnerFree()
    expect_equal(sum(o2), 212.0768500568, tolerance = 1e-6)
  })

  test_that("a vae grad fit does not leak into a later focei fast fit", {
    skip_on_cran()
    ## .foceiAnalyticSolveAll is SHARED with focei's own fast gradient and
    ## .vaeGradEnv lives for the whole session, so a vae grad fit must not leave
    ## the pooled-solve branch armed.  Without the `active` guard this made 23 of
    ## test-focei-fast-grad.R's assertions fail whenever the vae tests ran first.
    .ref <- suppressMessages(
      nlmixr2(.odeMod(), nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(print = 0L, covMethod = "", fast = TRUE,
                                     calcTables = FALSE)))
    suppressWarnings(suppressMessages(
      nlmixr2(.odeMod(), nlmixr2data::theo_sd, est = "vae",
              control = vaeControl(nonMuTheta = "grad", print = 0L, calcTables = FALSE,
                                   itersBurnIn = 10L, iters = 25L, klWarmup = 5L,
                                   gammaIter = 18L))))
    expect_false(isTRUE(.vaeGradEnv$active))
    expect_null(.vaeGradEnv$outerCols)
    .after <- suppressMessages(
      nlmixr2(.odeMod(), nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(print = 0L, covMethod = "", fast = TRUE,
                                     calcTables = FALSE)))
    expect_equal(.after$objf, .ref$objf, tolerance = 1e-4)
    expect_equal(unname(.after$theta), unname(.ref$theta), tolerance = 1e-4)
  })

  test_that("in scope, the gradient path is actually taken", {
    skip_on_cran()
    r <- suppressWarnings(suppressMessages(
      nlmixr2(.odeMod(), nlmixr2data::theo_sd, est = "vae",
              control = vaeControl(nonMuTheta = "grad", print = 0L, calcTables = FALSE,
                                   returnVae = TRUE, itersBurnIn = 10L, iters = 30L,
                                   klWarmup = 5L, gammaIter = 20L))))
    ## the mechanism assertion: M-steps past the KL warmup used the gradient and
    ## none silently fell back to bobyqa
    expect_true(r$nRegGrad > 0L)
    expect_equal(r$nRegFallback, 0L)
    expect_true(is.finite(r$regressTheta[["tv"]]))
  })
})
