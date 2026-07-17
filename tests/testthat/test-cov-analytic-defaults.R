# Analytic-first covMethod defaults across the mixed-model estimation methods.
# Always-run core unit test: no fits, only the control objects and the shared
# base-est mapping used by the post-fit covariance recompute.

test_that("covMethod defaults prefer analytic", {
  ## saem keeps the stochastic-approximation FIM ("sa") first, analytic second
  expect_identical(saemControl()$covMethod, "sa")
  expect_identical(saemControl(covMethod = "analytic")$covMethod, "analytic")
  expect_identical(nlmeControl()$covMethod, "analytic")
  expect_identical(vaeControl()$covMethod, "analytic")
  ## advi keeps its own variational covariance first, with analytic second
  expect_identical(adviControl()$covMethod, "advi")
  expect_identical(adviControl(covMethod = "analytic")$covMethod, "analytic")
  ## focei was already analytic-first
  expect_identical(foceiControl()$covMethod, 2L)
  expect_identical(foceiControl()$covType, "analytic")
  ## imp family inherits foceiControl's analytic default (integer slot + covType)
  expect_identical(impmapControl()$covMethod, 2L)
  expect_identical(impmapControl()$covType, "analytic")
})

test_that("nlme covMethod tokens round-trip", {
  expect_identical(nlmeControl(covMethod = "nlme")$covMethod, "nlme")
  expect_identical(nlmeControl(covMethod = "")$covMethod, "")
  expect_identical(nlmeControl(covMethod = "r,s")$covMethod, "r,s")
  expect_error(nlmeControl(covMethod = "bogus"))
})

test_that("saem keeps its other covMethod choices", {
  for (cm in c("linFim", "fim", "sa", "r,s", "r", "s")) {
    expect_identical(saemControl(covMethod = cm)$covMethod, cm)
  }
  expect_identical(saemControl(covMethod = "")$covMethod, "")
})

test_that(".foceiRecomputeBaseEst maps the covariance-recompute methods", {
  ## mu/irls families recompute on their own base method
  expect_identical(nlmixr2est:::.foceiRecomputeBaseEst("mfocei"), "focei")
  expect_identical(nlmixr2est:::.foceiRecomputeBaseEst("ifoce"), "foce")
  expect_identical(nlmixr2est:::.foceiRecomputeBaseEst("mfocep"), "focep")
  ## EM / nonparametric / nlme recompute on a zero-iteration focei model
  for (e in c("imp", "impmap", "qrpem", "nlme",
              "npag", "npb", "mnpag", "inpag", "mnpb", "inpb")) {
    expect_identical(nlmixr2est:::.foceiRecomputeBaseEst(e), "focei")
  }
  ## methods with their own covariance (or none) are not recomputed
  for (e in c("focei", "saem", "fsaem", "nlm", "nls")) {
    expect_null(nlmixr2est:::.foceiRecomputeBaseEst(e))
  }
})
