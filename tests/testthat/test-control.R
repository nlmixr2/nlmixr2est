nlmixrControlTest <- function(control) {
  # Solving options for table
  expect_true(inherits(control$rxControl, "rxControl"))
  # Options needed for parameter table generation
  expect_true(checkmate::testNumeric(control$ci, lower=0, upper=1, any.missing=FALSE, len=1))
  expect_true(checkmate::testIntegerish(control$sigdigTable, lower=1, any.missing=FALSE, len=1))
  expect_true(checkmate::testLogical(control$genRxControl, any.missing=FALSE, len=1))
  expect_true(checkmate::testLogical(control$calcTables, any.missing=FALSE, len=1))
  expect_true(checkmate::testLogical(control$compress, any.missing=FALSE, len=1))
}

test_that("test foceiControl option sanity", {
  expect_error(foceiControl(), NA)
  nlmixrControlTest(foceiControl())

  .ctl <- foceiControl()
  expect_error(do.call(foceiControl, .ctl), NA)
  .ctl2 <- do.call(foceiControl, .ctl)
  expect_equal(.ctl, .ctl2)

  # ResetEtaP
  .ctl <- foceiControl(resetEtaP=0.5)
  .ctl2 <- do.call(foceiControl, .ctl)
  expect_equal(.ctl, .ctl2)

  # resetThetaP
  .ctl <- foceiControl(resetThetaP=0.5)
  .ctl2 <- do.call(foceiControl, .ctl)
  expect_equal(.ctl, .ctl2)

  # resetThetaFinalP
  .ctl <- foceiControl(resetThetaFinalP=0.5)
  .ctl2 <- do.call(foceiControl, .ctl)
  expect_equal(.ctl, .ctl2)
  expect_true(.ctl$genRxControl)

  .ctl <- foceiControl(rxControl=rxControl(sigdig=6))
  expect_false(.ctl$genRxControl)
  .ctl2 <- do.call(foceiControl, .ctl)
  expect_equal(.ctl, .ctl2)

  expect_error(foceiControl(foceiControl="matt"))
})


test_that("saemControl sanity", {
  expect_error(saemControl(), NA)
  nlmixrControlTest(saemControl())
  .ctl <- saemControl()
  expect_error(do.call(saemControl, .ctl), NA)
  .ctl2 <- do.call(saemControl, .ctl)
  expect_equal(.ctl, .ctl2)

  .ctl <- saemControl(rxControl=rxControl(sigdig=6))
  expect_false(.ctl$genRxControl)
  .ctl2 <- do.call(saemControl, .ctl)
  expect_equal(.ctl, .ctl2)

  .ctl <- saemControl(trace=1)
  .ctl2 <- do.call(saemControl, .ctl)
  expect_equal(.ctl, .ctl2)

  expect_error(saemControl(foceiControl="matt"))
})

test_that("nlmixr2NlmeControl sanity", {
  expect_error(nlmixr2NlmeControl(), NA)
  nlmixrControlTest(nlmixr2NlmeControl())

  .ctl <- nlmixr2NlmeControl()
  expect_error(do.call(nlmixr2NlmeControl, .ctl), NA)
  .ctl2 <- do.call(nlmixr2NlmeControl, .ctl)
  expect_equal(.ctl, .ctl2)

  expect_error(nlmixr2NlmeControl(foceiControl="matt"))
})

test_that("foceiControl for lbfgsb3c", {
  .tmp <- foceiControl(print = 1, outerOpt="lbfgsb3c")
  expect_error(do.call("foceiControl", .tmp), NA)
  .tmp2 <- do.call("foceiControl", .tmp)
  expect_equal(.tmp, .tmp2)
})

test_that("saemControl can take integer for covMethod", {
  expect_error(saemControl(covMethod=0L), NA)
})

test_that("saemControl can take integer for covMethod", {
  expect_error(nlmeControl(covMethod=0L), NA)
})
