## #1051: foceiControl(fast=TRUE)'s analytic outer gradient re-derives
## d(objective)/d(theta) from MODEL sensitivities only, so it cannot carry an
## external likelihood contribution's own theta dependence.  It must decline to
## finite differences for a contributor that CHANGES the objective, and stay on
## the fast path for a pure observer (which changes nothing and is exact).
##
## The stock test contributor's constant per-observation shift cannot move an
## optimum, so it cannot show this; `_nlmixr2est_setTestContribAddLLf` adds c*f,
## which reaches theta through the prediction and supplies its own exact
## d(LL)/d(eta) -- so the inner mode stays exact and any outer difference is
## attributable to the outer gradient alone.

test_that("fast=TRUE declines the analytic gradient for a theta-dependent contributor", {
  skip_on_cran()

  .old <- rxode2::getRxThreads()
  on.exit(rxode2::setRxThreads(.old), add = TRUE)
  rxode2::setRxThreads(1L)   # test contributor uses global accumulators

  d <- data.frame(ID = rep(1:3, each = 3),
                  TIME = rep(c(0.5, 1, 2), 3),
                  DV = c(1.4, 1.5, 1.3, 1.7, 1.2, 1.6, 1.1, 1.9, 1.5),
                  EVID = 0L)

  mod <- function() {
    ini({
      level <- 0.3
      etaLevel ~ 0.1
      error <- fix(0.2)
    })
    model({
      prediction <- exp(level + etaLevel)
      prediction ~ prop(error)
    })
  }

  .run <- function(fast, cc) {
    .Call("_nlmixr2est_registerTestContrib", PACKAGE = "nlmixr2est")
    on.exit(.Call("_nlmixr2est_removeTestContrib", PACKAGE = "nlmixr2est"), add = TRUE)
    .Call("_nlmixr2est_setTestContribAddLLf", cc, PACKAGE = "nlmixr2est")
    suppressWarnings(.nlmixr(mod, d, est = "focei",
                             control = foceiControl(print = 0L, fast = fast,
                                                    calcTables = FALSE)))
  }

  ## observer only (c = 0): the contributor writes nothing back, so the analytic
  ## gradient is exact and must still be used
  .obsT <- .run(TRUE, 0)
  .obsF <- .run(FALSE, 0)
  expect_gt(as.integer(.obsT$env$nAnalyticGradDirect), 0L)
  expect_equal(as.integer(.obsT$env$nFDGradFast), 0L)
  expect_equal(unname(fixef(.obsT)["level"]), unname(fixef(.obsF)["level"]),
               tolerance = 1e-3)

  ## theta-dependent contribution: the analytic gradient must decline
  .cT <- .run(TRUE, -0.5)
  .cF <- .run(FALSE, -0.5)
  expect_equal(as.integer(.cT$env$nAnalyticGradDirect), 0L)
  expect_gt(as.integer(.cT$env$nFDGradFast), 0L)

  ## ...and having declined, fast= now selects the gradient and nothing else:
  ## before the fix the two arms differed by 5.5e-02 in `level` and 0.88 in objf
  expect_equal(unname(fixef(.cT)["level"]), unname(fixef(.cF)["level"]),
               tolerance = 1e-4)
  expect_equal(.cT$objf, .cF$objf, tolerance = 1e-5)

  ## the contribution really did move the optimum (so the comparison is not vacuous)
  expect_gt(abs(unname(fixef(.cF)["level"]) - unname(fixef(.obsF)["level"])), 1e-3)
})
