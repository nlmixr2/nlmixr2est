nmTest({
  nlmixrControlTest <- function(control) {
    # Solving options for table
    expect_true(inherits(control$rxControl, "rxControl"))
    # Options needed for parameter table generation
    expect_true(checkmate::testNumeric(control$ci, lower = 0, upper = 1, any.missing = FALSE, len = 1))
    expect_true(checkmate::testIntegerish(control$sigdigTable, lower = 1, any.missing = FALSE, len = 1))
    expect_true(checkmate::testLogical(control$genRxControl, any.missing = FALSE, len = 1))
    expect_true(checkmate::testLogical(control$calcTables, any.missing = FALSE, len = 1))
    expect_true(checkmate::testLogical(control$compress, any.missing = FALSE, len = 1))
  }

  test_that("test foceiControl option sanity", {
    expect_error(foceiControl(), NA)
    nlmixrControlTest(foceiControl())

    .ctl <- foceiControl()
    expect_error(do.call(foceiControl, .ctl), NA)
    .ctl2 <- do.call(foceiControl, .ctl)
    expect_equal(.ctl, .ctl2)

    # ResetEtaP
    .ctl <- foceiControl(resetEtaP = 0.5)
    .ctl2 <- do.call(foceiControl, .ctl)
    expect_equal(.ctl, .ctl2)

    # resetThetaP
    .ctl <- foceiControl(resetThetaP = 0.5)
    .ctl2 <- do.call(foceiControl, .ctl)
    expect_equal(.ctl, .ctl2)

    # resetThetaFinalP
    .ctl <- foceiControl(resetThetaFinalP = 0.5)
    .ctl2 <- do.call(foceiControl, .ctl)
    expect_equal(.ctl, .ctl2)
    expect_true(.ctl$genRxControl)

    .ctl <- foceiControl(rxControl = rxControl(sigdig = 6))
    expect_false(.ctl$genRxControl)
    .ctl2 <- do.call(foceiControl, .ctl)
    expect_equal(.ctl, .ctl2)

    expect_error(foceiControl(foceiControl = "matt"))
  })

  test_that("saemControl sanity", {
    expect_error(saemControl(), NA)
    nlmixrControlTest(saemControl())
    .ctl <- saemControl()
    expect_error(do.call(saemControl, .ctl), NA)
    .ctl2 <- do.call(saemControl, .ctl)
    expect_equal(.ctl, .ctl2)

    .ctl <- saemControl(rxControl = rxControl(sigdig = 6))
    expect_false(.ctl$genRxControl)
    .ctl2 <- do.call(saemControl, .ctl)
    expect_equal(.ctl, .ctl2)

    .ctl <- saemControl(trace = 1)
    .ctl2 <- do.call(saemControl, .ctl)
    expect_equal(.ctl, .ctl2)

    ## covMethod="" requests no covariance; it is a documented choice but
    ## match.arg() cannot select it (pmatch("") matches nothing), so it must be
    ## handled explicitly rather than erroring.
    expect_error(saemControl(covMethod = ""), NA)
    expect_equal(saemControl(covMethod = "")$covMethod, "")
    .ctl <- saemControl(covMethod = "")
    .ctl2 <- do.call(saemControl, .ctl)
    expect_equal(.ctl, .ctl2)

    expect_error(saemControl(foceiControl = "matt"))
  })

  test_that("saemControl rejects the withdrawn lbfgs* options (#878)", {
    ## announced in 7.0.2 but never implemented; they must not silently be
    ## accepted and stored again
    expect_error(saemControl(lbfgsLmm = 5L), "unused argument")
    expect_error(saemControl(lbfgsFactr = 1e7), "unused argument")
    expect_error(saemControl(lbfgsPgtol = 0), "unused argument")
    expect_error(saemControl(lbfgsMaxIter = 20L), "unused argument")
    expect_false(any(c("lbfgsLmm", "lbfgsFactr", "lbfgsPgtol", "lbfgsMaxIter") %in% names(saemControl())))
  })

  test_that("nlmixr2NlmeControl sanity", {
    expect_error(nlmixr2NlmeControl(), NA)
    nlmixrControlTest(nlmixr2NlmeControl())

    .ctl <- nlmixr2NlmeControl()
    expect_error(do.call(nlmixr2NlmeControl, .ctl), NA)
    .ctl2 <- do.call(nlmixr2NlmeControl, .ctl)
    expect_equal(.ctl, .ctl2)

    expect_error(nlmixr2NlmeControl(foceiControl = "matt"))
  })

  test_that("foceiControl for lbfgsb3c", {
    .tmp <- foceiControl(print = 1, outerOpt = "lbfgsb3c")
    expect_error(do.call("foceiControl", .tmp), NA)
    .tmp2 <- do.call("foceiControl", .tmp)
    expect_equal(.tmp, .tmp2)
  })

  test_that("saemControl can take integer for covMethod", {
    expect_error(saemControl(covMethod = 0L), NA)
  })

  test_that("saemControl can take integer for covMethod", {
    expect_error(nlmeControl(covMethod = 0L), NA)
  })
})

nmTest({
  test_that("saemControl(nu=) survives the plain-list rebuild saem does", {
    ## .saemFamilyControl() rebuilds the control with
    ## do.call(saemControl, <plain list>), which routes `mcmc` through the
    ## `.xtra$mcmc` branch.  That branch copied nBurn, nEm and nmc out of it but
    ## only VALIDATED nu, so every saemControl(nu = ...) silently fitted with
    ## the default c(2, 2, 2).
    .c <- saemControl(nBurn = 7, nEm = 9, nmc = 4, nu = c(10, 11, 12))
    .r <- do.call(saemControl, unclass(.c))
    expect_equal(.r$mcmc$nu, c(10, 11, 12))
    expect_equal(.r$mcmc$niter, c(7, 9))
    expect_equal(.r$mcmc$nmc, 4)
  })
})

test_that("saemControl() checks nSaCov and seed as saControl() does", {
  .msg <- "Assertion on 'nSaCov' failed: Element 1 is not >= 1."
  expect_error(saemControl(nSaCov = 0), .msg, fixed = TRUE)
  expect_error(saControl(nSaCov = 0), .msg, fixed = TRUE)
  expect_error(saemControl(nSaCov = -5), .msg, fixed = TRUE)
  expect_error(
    saemControl(nSaCov = NA),
    "Assertion on 'nSaCov' failed: Contains missing values (element 1).",
    fixed = TRUE
  )
  expect_error(
    saemControl(nSaCov = 1.5),
    "Assertion on 'nSaCov' failed: Must be of type 'integerish', but element 1 is not close to an integer.",
    fixed = TRUE
  )
  expect_error(
    saemControl(nSaCov = c(10, 20)),
    "Assertion on 'nSaCov' failed: Must have length 1, but has length 2.",
    fixed = TRUE
  )
  expect_identical(saemControl(nSaCov = 1)$nSaCov, 1L)
  # one seed: the fit seeds its sampler from a single value
  .msg <- "Assertion on 'seed' failed: Must have length 1, but has length 2."
  expect_error(saemControl(seed = c(1, 2)), .msg, fixed = TRUE)
  expect_error(saControl(seed = c(1, 2)), .msg, fixed = TRUE)
  expect_identical(saemControl(seed = 3)$seed, 3)
})
