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
