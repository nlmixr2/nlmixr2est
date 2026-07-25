# viControl() argument handling + est="emvi"/"fbvi" dispatch resolution.
# Always-run core unit test: no fit, just the control object and that dispatch
# reaches the two variational methods.

test_that("viControl() builds and validates", {
  ctl <- viControl()
  expect_s3_class(ctl, "viControl")
  expect_identical(ctl$viFamily, "fullRank")
  ## pointEstimate is resolved from `est`, not defaulted in the control
  expect_null(ctl$pointEstimate)
  expect_identical(ctl$optim, "advi")
  expect_identical(ctl$seed, 42L)

  ## enum + range validation
  expect_error(viControl(viFamily = "bogus"))
  expect_error(viControl(iters = 0L))
  expect_error(viControl(alpha = 2))          # alpha in (0,1)
  expect_error(viControl(nMc = 0L))
  expect_error(viControl(bogusArg = 1))       # unused argument

  ## knobs round-trip
  ctl2 <- viControl(iters = 50L, nMc = 3L, viFamily = "meanField",
                      pointEstimate = FALSE, optim = "adam")
  expect_identical(ctl2$iters, 50L)
  expect_identical(ctl2$nMc, 3L)
  expect_identical(ctl2$viFamily, "meanField")
  expect_false(ctl2$pointEstimate)
  expect_identical(ctl2$optim, "adam")
})

test_that("getValidNlmixrCtl.emvi/.fbvi normalize control and resolve pointEstimate", {
  v <- getValidNlmixrCtl.emvi(list(viControl(iters = 7L)))
  expect_s3_class(v, "viControl")
  expect_identical(v$iters, 7L)
  ## a bare list is coerced to viControl
  v2 <- getValidNlmixrCtl.emvi(list(list(iters = 9L)))
  expect_s3_class(v2, "viControl")
  expect_identical(v2$iters, 9L)

  ## the two methods differ only in the inference axis, resolved from `est`
  expect_true(getValidNlmixrCtl.emvi(list(viControl()))$pointEstimate)
  expect_false(getValidNlmixrCtl.fbvi(list(viControl()))$pointEstimate)

  ## an explicit pointEstimate that AGREES with est is accepted
  expect_true(getValidNlmixrCtl.emvi(list(viControl(pointEstimate = TRUE)))$pointEstimate)
  expect_false(getValidNlmixrCtl.fbvi(list(viControl(pointEstimate = FALSE)))$pointEstimate)

  ## one that CONTRADICTS it loses to `est` -- but is ANNOUNCED, not silent.
  ## `est` has to win rather than error because the two methods share one control
  ## class, so re-estimating a fit with the other method pipes the completed
  ## fit's control (carrying the old pointEstimate) forward.
  expect_message(v3 <- getValidNlmixrCtl.emvi(list(viControl(pointEstimate = FALSE))),
                 "pointEstimate=TRUE")
  expect_true(v3$pointEstimate)
  expect_message(v4 <- getValidNlmixrCtl.fbvi(list(viControl(pointEstimate = TRUE))),
                 "pointEstimate=FALSE")
  expect_false(v4$pointEstimate)
  ## an agreeing value says nothing
  expect_silent(getValidNlmixrCtl.emvi(list(viControl(pointEstimate = TRUE))))
})

test_that("rxUiDeparse.viControl round-trips changed args", {
  dp <- rxUiDeparse(viControl(iters = 50L, viFamily = "meanField"), "x")
  expect_true(any(grepl("viControl", as.character(dp))))
  expect_true(any(grepl("meanField", as.character(dp))))
})

test_that("perNoCor takes a fraction or an absolute iteration count", {
  ## <= 1 is a fraction of the run
  expect_equal(viControl(perNoCor = 0.5)$perNoCor, 0.5)
  expect_equal(viControl(perNoCor = 1)$perNoCor, 1)
  ## > 1 is an ABSOLUTE iteration count -- the only form a resumed fit can
  ## reproduce, so it has to survive the control (it used to be capped at 1)
  expect_equal(viControl(perNoCor = 90)$perNoCor, 90)
  expect_equal(vaeControl(perNoCor = 90)$perNoCor, 90)
  ## a fractional absolute count is rejected rather than silently rounded to a
  ## schedule the user did not ask for
  expect_error(viControl(perNoCor = 2.5), "perNoCor")
  expect_error(vaeControl(perNoCor = 2.5), "perNoCor")
  expect_error(viControl(perNoCor = -1), "perNoCor")
})
