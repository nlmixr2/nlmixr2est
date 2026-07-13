# adviControl() argument handling + est="advi" dispatch resolution.  Always-run
# core unit test: no fit, just the control object and that dispatch reaches the
# advi method.

test_that("adviControl() builds and validates", {
  ctl <- adviControl()
  expect_s3_class(ctl, "adviControl")
  expect_identical(ctl$adviFamily, "fullRank")
  expect_true(ctl$pointEstimate)
  expect_identical(ctl$optim, "advi")
  expect_identical(ctl$seed, 42L)

  ## enum + range validation
  expect_error(adviControl(adviFamily = "bogus"))
  expect_error(adviControl(iters = 0L))
  expect_error(adviControl(alpha = 2))          # alpha in (0,1)
  expect_error(adviControl(nMc = 0L))
  expect_error(adviControl(bogusArg = 1))       # unused argument

  ## knobs round-trip
  ctl2 <- adviControl(iters = 50L, nMc = 3L, adviFamily = "meanField",
                      pointEstimate = FALSE, optim = "adam")
  expect_identical(ctl2$iters, 50L)
  expect_identical(ctl2$nMc, 3L)
  expect_identical(ctl2$adviFamily, "meanField")
  expect_false(ctl2$pointEstimate)
  expect_identical(ctl2$optim, "adam")
})

test_that("getValidNlmixrCtl.advi normalizes control", {
  v <- getValidNlmixrCtl.advi(list(adviControl(iters = 7L)))
  expect_s3_class(v, "adviControl")
  expect_identical(v$iters, 7L)
  ## a bare list is coerced to adviControl
  v2 <- getValidNlmixrCtl.advi(list(list(iters = 9L)))
  expect_s3_class(v2, "adviControl")
  expect_identical(v2$iters, 9L)
})

test_that("rxUiDeparse.adviControl round-trips changed args", {
  dp <- rxUiDeparse(adviControl(iters = 50L, adviFamily = "meanField"), "x")
  expect_true(any(grepl("adviControl", as.character(dp))))
  expect_true(any(grepl("meanField", as.character(dp))))
})
