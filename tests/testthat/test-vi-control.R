# emviControl()/fbviControl() argument handling + est="emvi"/"fbvi" dispatch
# resolution.  Always-run core unit test: no fit, just the control objects and
# that dispatch reaches the two variational methods.

test_that("emviControl() builds and validates", {
  ctl <- emviControl()
  expect_s3_class(ctl, "emviControl")
  expect_identical(ctl$viFamily, "fullRank")
  ## pointEstimate is resolved from `est`, not defaulted in the control
  expect_null(ctl$pointEstimate)
  expect_identical(ctl$optim, "advi")
  expect_identical(ctl$seed, 42L)

  ## enum + range validation
  expect_error(emviControl(viFamily = "bogus"))
  expect_error(emviControl(iters = 0L))
  expect_error(emviControl(alpha = 2))          # alpha in (0,1)
  expect_error(emviControl(nMc = 0L))
  expect_error(emviControl(bogusArg = 1))       # unused argument

  ## knobs round-trip
  ctl2 <- emviControl(iters = 50L, nMc = 3L, viFamily = "meanField",
                      pointEstimate = FALSE, optim = "adam")
  expect_identical(ctl2$iters, 50L)
  expect_identical(ctl2$nMc, 3L)
  expect_identical(ctl2$viFamily, "meanField")
  expect_false(ctl2$pointEstimate)
  expect_identical(ctl2$optim, "adam")
})

test_that("getValidNlmixrCtl.emvi/.fbvi normalize control and resolve pointEstimate", {
  v <- getValidNlmixrCtl.emvi(list(emviControl(iters = 7L)))
  expect_s3_class(v, "emviControl")
  expect_identical(v$iters, 7L)
  ## a bare list is coerced to emviControl
  v2 <- getValidNlmixrCtl.emvi(list(list(iters = 9L)))
  expect_s3_class(v2, "emviControl")
  expect_identical(v2$iters, 9L)

  ## the two methods differ only in the inference axis, resolved from `est`
  expect_true(getValidNlmixrCtl.emvi(list(emviControl()))$pointEstimate)
  expect_false(getValidNlmixrCtl.fbvi(list(emviControl()))$pointEstimate)

  ## an explicit pointEstimate that AGREES with est is accepted
  expect_true(getValidNlmixrCtl.emvi(list(emviControl(pointEstimate = TRUE)))$pointEstimate)
  expect_false(getValidNlmixrCtl.fbvi(list(emviControl(pointEstimate = FALSE)))$pointEstimate)

  ## one that CONTRADICTS it loses to `est` -- but is ANNOUNCED, not silent.
  ## `est` has to win rather than error because the two methods share one control
  ## class, so re-estimating a fit with the other method pipes the completed
  ## fit's control (carrying the old pointEstimate) forward.
  expect_message(v3 <- getValidNlmixrCtl.emvi(list(emviControl(pointEstimate = FALSE))),
                 "pointEstimate=TRUE")
  expect_true(v3$pointEstimate)
  expect_message(v4 <- getValidNlmixrCtl.fbvi(list(emviControl(pointEstimate = TRUE))),
                 "pointEstimate=FALSE")
  expect_false(v4$pointEstimate)
  ## an agreeing value says nothing
  expect_silent(getValidNlmixrCtl.emvi(list(emviControl(pointEstimate = TRUE))))
})

test_that("rxUiDeparse.emviControl round-trips changed args", {
  dp <- rxUiDeparse(emviControl(iters = 50L, viFamily = "meanField"), "x")
  expect_true(any(grepl("emviControl", as.character(dp))))
  expect_true(any(grepl("meanField", as.character(dp))))
})

test_that("perNoCor takes a fraction or an absolute iteration count", {
  ## <= 1 is a fraction of the run
  expect_equal(emviControl(perNoCor = 0.5)$perNoCor, 0.5)
  expect_equal(emviControl(perNoCor = 1)$perNoCor, 1)
  ## > 1 is an ABSOLUTE iteration count -- the only form a resumed fit can
  ## reproduce, so it has to survive the control (it used to be capped at 1)
  expect_equal(emviControl(perNoCor = 90)$perNoCor, 90)
  expect_equal(vaeControl(perNoCor = 90)$perNoCor, 90)
  ## a fractional absolute count is rejected rather than silently rounded to a
  ## schedule the user did not ask for
  expect_error(emviControl(perNoCor = 2.5), "perNoCor")
  expect_error(vaeControl(perNoCor = 2.5), "perNoCor")
  expect_error(emviControl(perNoCor = -1), "perNoCor")
})

test_that("fbviControl() is emviControl() with pointEstimate = FALSE", {
  ## same shim shape as impControl()/qrpemControl() over impmapControl(): a thin
  ## wrapper that sets the distinguishing field and keeps the base class, so one
  ## control function exists per est= value
  ctl <- fbviControl()
  expect_s3_class(ctl, "emviControl")
  expect_false(ctl$pointEstimate)

  ## every other argument still reaches emviControl()
  ctl2 <- fbviControl(iters = 33L, viFamily = "meanField")
  expect_identical(ctl2$iters, 33L)
  expect_identical(ctl2$viFamily, "meanField")
  expect_false(ctl2$pointEstimate)
  ## and its validation is emviControl()'s
  expect_error(fbviControl(viFamily = "bogus"))
  expect_error(fbviControl(bogusArg = 1))

  ## pointEstimate is a NAMED formal, so passing it explicitly binds rather than
  ## colliding as a duplicate argument (the qrpemControl(qr=) shape)
  expect_true(fbviControl(pointEstimate = TRUE)$pointEstimate)
  expect_error(fbviControl(pointEstimate = "yes"))

  ## POSITIONAL arguments must land on the same formals they would in
  ## emviControl().  Passing pointEstimate= into emviControl() from the shim
  ## would bind formal #5 by name, dropping it out of positional matching and
  ## shifting every later positional argument one slot -- "adam" would arrive at
  ## adaptEta.  seed/iters/nMc/viFamily/pointEstimate/optim are formals 1-6.
  p <- fbviControl(42L, 300L, 1L, "fullRank", FALSE, "adam")
  expect_identical(p$seed, 42L)
  expect_identical(p$iters, 300L)
  expect_identical(p$nMc, 1L)
  expect_identical(p$viFamily, "fullRank")
  expect_identical(p$optim, "adam")
  expect_true(p$adaptEta)                 # NOT "adam"
  expect_false(p$pointEstimate)

  ## defaults differ ONLY on that axis
  .drop <- function(x) { x$pointEstimate <- NULL; x }
  expect_equal(.drop(fbviControl()), .drop(emviControl()))
})

test_that("getValidNlmixrCtl.fbvi accepts NULL and a bare list", {
  expect_false(getValidNlmixrCtl.fbvi(list(NULL))$pointEstimate)
  v <- getValidNlmixrCtl.fbvi(list(list(iters = 9L)))
  expect_s3_class(v, "emviControl")
  expect_identical(v$iters, 9L)
  expect_false(v$pointEstimate)
  ## an fbviControl() passed to est="emvi" is overridden by est, with a message
  expect_message(v2 <- getValidNlmixrCtl.emvi(list(fbviControl())), "pointEstimate=TRUE")
  expect_true(v2$pointEstimate)
})
