# Per-family covMethod defaults across the mixed-model estimation methods.
# Always-run core unit test: no fits, only the control objects and the shared
# base-est mapping used by the post-fit covariance recompute.

test_that("per-family covMethod defaults", {
  ## saem keeps the stochastic-approximation FIM ("sa")
  expect_identical(saemControl()$covMethod, "sa")
  expect_identical(saemControl(covMethod = "analytic")$covMethod, "analytic")
  ## nlme now keeps nlme's own covariance by default
  expect_identical(nlmeControl()$covMethod, "nlme")
  ## vae now defaults to the r,s sandwich
  expect_identical(vaeControl()$covMethod, "r,s")
  ## advi keeps its own variational covariance first, with analytic second
  expect_identical(adviControl()$covMethod, "advi")
  expect_identical(adviControl(covMethod = "analytic")$covMethod, "analytic")
  ## focei now defaults to the r,s sandwich (integer slot 1, finite-difference)
  expect_identical(foceiControl()$covMethod, 1L)
  expect_identical(foceiControl()$covType, "fd")
  expect_identical(foceiControl(covMethod = "analytic")$covMethod, 2L)
  ## imp family keeps the importance-sampling default (maps to the analytic slot
  ## internally, with the MC covariance driven by the impCov flag)
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
