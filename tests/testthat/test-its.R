# est="its" controls (no model fits); fit tests are in test-its-fit.R
nmTest({
  test_that("itsControl() sets the deterministic ITS kernel", {
    .ctl <- itsControl()
    expect_s3_class(.ctl, "impmapControl")
    expect_true(.ctl$its)
    expect_identical(.ctl$nIter, 50L)
    expect_identical(.ctl$ctol, 1e-4)
    expect_identical(.ctl$isample, 1L)
    expect_identical(.ctl$mapIter, 1L)
    expect_identical(.ctl$nConvWindow, 1L)
    expect_false(.ctl$auto)
    expect_false(.ctl$qr)
    expect_false(.ctl$sir)
    expect_identical(.ctl$iaccept, 0)
    expect_identical(.ctl$interaction, 1L)
    expect_identical(itsControl(interaction = FALSE)$interaction, 0L)
    expect_identical(itsControl(maxIter = 7, tol = 1e-3)$nIter, 7L)
    expect_identical(itsControl(maxIter = 7, tol = 1e-3)$ctol, 1e-3)
    # foceiControl options pass through
    expect_identical(itsControl(maxInnerIterations = 20L)$maxInnerIterations, 20L)
    expect_error(itsControl(maxIter = -1))
    expect_error(itsControl(tol = -1))
    expect_error(itsControl(covMethod = "imp"))
  })

  test_that("itsControl() round-trips and cannot be talked out of ITS", {
    .ctl <- itsControl(maxIter = 12L, tol = 1e-5, interaction = FALSE, covMethod = "r")
    .v <- getValidNlmixrCtl.its(structure(list(.ctl), class = "its"))
    expect_true(.v$its)
    expect_identical(.v$nIter, 12L)
    expect_identical(.v$ctol, 1e-5)
    expect_identical(.v$interaction, 0L)
    expect_identical(.v$covMethod, .ctl$covMethod)
    # a sampling control from another method keeps only the shared options
    .v <- getValidNlmixrCtl.its(structure(list(impmapControl(isample = 500L, maxInnerIterations = 30L)), class = "its"))
    expect_true(.v$its)
    expect_identical(.v$isample, 1L)
    expect_identical(.v$nIter, 50L)
    expect_identical(.v$maxInnerIterations, 30L)
    .v <- getValidNlmixrCtl.its(structure(list(NULL), class = "its"))
    expect_true(.v$its)
    .v <- getValidNlmixrCtl.its(structure(list(list(maxIter = 3L)), class = "its"))
    expect_identical(.v$nIter, 3L)
  })

  test_that("an ITS control used by impmap gets the sampling defaults back", {
    .def <- impmapControl()
    .v <- suppressMessages(getValidNlmixrCtl.impmap(structure(list(itsControl()), class = "impmap")))
    expect_null(.v$its)
    expect_identical(.v$isample, .def$isample)
    expect_identical(.v$nIter, .def$nIter)
    expect_identical(.v$iaccept, .def$iaccept)
    expect_identical(.v$auto, .def$auto)
    expect_message(
      getValidNlmixrCtl.qrpem(structure(list(itsControl()), class = "qrpem")),
      "importance-sampling defaults"
    )
  })

  test_that("an ITS control used by qrpem or the np engines is converted", {
    .v <- suppressMessages(getValidNlmixrCtl.qrpem(structure(list(itsControl()), class = "qrpem")))
    expect_true(.v$qr)
    expect_true(.v$sir)
    expect_identical(.v$isample, qrpemControl()$isample)
    for (.e in c("npag", "npb")) {
      .v <- suppressMessages(getValidNlmixrControl(itsControl(), .e))
      expect_null(.v$its)
      expect_identical(.v$est, .e)
      .ctl <- if (.e == "npag") npagControl(cycles = 3L) else npbControl()
      .v <- suppressMessages(getValidNlmixrControl(.ctl, "its"))
      expect_true(.v$its)
      expect_null(.v$cycles)
    }
  })

  test_that("itsControl() deparses to itself", {
    expect_false(grepl("sirSample|nConvWindow|ctol", deparse1(rxUiDeparse(itsControl(), "ctl"))))
    .ctl <- itsControl(maxIter = 7L, tol = 1e-3, interaction = FALSE, maxInnerIterations = 20L)
    .d <- rxUiDeparse(.ctl, "ctl")
    expect_match(deparse1(.d), "itsControl")
    .ctl2 <- eval(.d[[3]])
    expect_true(.ctl2$its)
    expect_identical(.ctl2$nIter, 7L)
    expect_identical(.ctl2$ctol, 1e-3)
    expect_identical(.ctl2$interaction, 0L)
    expect_identical(.ctl2$maxInnerIterations, 20L)
  })

  test_that("its is a registered estimation method", {
    expect_true("its" %in% nlmixr2AllEst())
    expect_true(is.function(getS3method("nlmixr2Est", "its")))
    expect_identical(.foceiRecomputeBaseEst("its"), "focei")
  })
})
