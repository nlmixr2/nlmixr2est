test_that("impmapControl option sanity and IS defaults", {
  expect_error(impmapControl(), NA)
  .ctl <- impmapControl()
  expect_s3_class(.ctl, "impmapControl")
  # inherits the FOCEI/MAP option surface
  expect_false(is.null(.ctl$maxOuterIterations))
  # mu-referencing is forced on for the MAP proposal center
  expect_identical(.ctl$muModel, "lin")
  # importance-sampling / EM defaults
  expect_identical(.ctl$isample, 300L)
  expect_identical(.ctl$nIter, 100L)
  expect_identical(.ctl$mapIter, 1L)
  expect_identical(.ctl$gamma, 1.0)
  expect_identical(.ctl$iscaleMin, 0.1)
  expect_identical(.ctl$iscaleMax, 10.0)
  expect_identical(.ctl$iaccept, 0.4)
  expect_null(.ctl$ctol)
  expect_identical(.ctl$nConvWindow, 10L)
  expect_identical(.ctl$impSeed, 42L)

  # round-trips through do.call (used by getValidNlmixrCtl / .foceiFamilyControl)
  expect_error(do.call(impmapControl, .ctl), NA)
})

test_that("impmapControl overrides and focei passthrough", {
  .ctl <- impmapControl(isample = 50L, gamma = 2, impSeed = 7L,
                        maxOuterIterations = 3L)
  expect_identical(.ctl$isample, 50L)
  expect_identical(.ctl$gamma, 2.0)
  expect_identical(.ctl$impSeed, 7L)
  # forwarded to foceiControl via ...
  expect_identical(.ctl$maxOuterIterations, 3L)
})

test_that("down-conversion to foceiControl strips IS-only names", {
  .env <- new.env()
  .env$impmapControl <- impmapControl()
  .fc <- nlmixr2est:::.impmapControlToFoceiControl(.env, assign = FALSE)
  expect_s3_class(.fc, "foceiControl")
  # IS/EM-only names must not leak into the plain foceiControl
  for (.n in nlmixr2est:::.impmapIsControlNames) {
    expect_null(.fc[[.n]])
  }
  # MAP-relevant focei options are preserved
  expect_identical(.fc$muModel, "lin")
})

test_that("impmap dispatch is registered and discoverable", {
  expect_true("impmap" %in% nlmixr2AllEst())
  expect_true(is.function(getS3method("nlmixr2Est", "impmap")))
  # mu-hook activation gate is a control-dependent predicate (mufocei-style)
  expect_true(is.function(attr(nlmixr2Est.impmap, "mu")))
})

test_that("getValidNlmixrCtl.impmap yields a default impmapControl", {
  expect_s3_class(getValidNlmixrCtl.impmap(list(NULL)), "impmapControl")
})

test_that("M0: impmap fit stops with an explicit under-construction error", {
  one.cmt <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  expect_error(
    nlmixr2(one.cmt, nlmixr2data::theo_sd, "impmap", impmapControl(print = 0L)),
    "under construction"
  )
})
