# Fit-based QRPEM tests -- registered in .slowBatches (tests/testthat.R), so
# they run in the weekly slow-tests workflow, not on every push/PR.
# Quick no-fit QRPEM tests live in test-qrpem.R.
nmTest({
  .oneCmt <- function() {
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

  test_that("Q1: qr/sir options round-trip into the C++ kernel", {
    .f <- suppressWarnings(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, "impmap",
              impmapControl(print=0L, nIter=1L, isample=50L,
                            sir=TRUE, sirSample=25L)))
    .e <- .f$env
    # qr not requested; sir requested: values must come back from op_focei
    expect_false(.e$impQr)
    expect_true(.e$impSir)
    expect_identical(.e$impSirSample, 25L)
  })

  test_that("Q1: default (qr/sir off) fit is unchanged vs the pre-QRPEM baseline", {
    .ref <- readRDS(test_path("fixtures", "qrpem-baseline-ref.rds"))
    .f <- suppressWarnings(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, "impmap",
              impmapControl(print=0L, nIter=5L, isample=100L)))
    expect_equal(fixef(.f), .ref$fixef, tolerance=1e-8)
    expect_equal(.f$omega, .ref$omega, tolerance=1e-8)
    expect_equal(.f$env$impObj, .ref$obj, tolerance=1e-8)
    # the E-step samples themselves are seed-pinned -> same draw stream
    expect_equal(.f$env$impSamples[[1]], .ref$samples1, tolerance=1e-8)
  })
})
