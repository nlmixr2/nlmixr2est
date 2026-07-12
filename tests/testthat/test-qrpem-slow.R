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

  # Back-transform a subject's stashed samples to the underlying N(0,1) points:
  # S = mode + L Z with L = lower chol(gamma * H^-1), so Z = L^-1 (S - mode).
  .zPoints <- function(env, i) {
    .g <- env$impGammaTrace[length(env$impGammaTrace)]
    .L <- t(chol(.g * solve(env$impEtaHess[[i]])))
    t(solve(.L, t(env$impSamples[[i]]) - as.numeric(env$impEtaMode[i, ])))
  }

  test_that("Q3: qr proposal matches N(mode, gamma*H^-1) and reweights to H^-1", {
    .gamma <- 2
    .f <- suppressWarnings(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, "impmap",
              impmapControl(print=0L, nIter=1L, isample=2048L, gamma=.gamma,
                            qr=TRUE)))
    .e <- .f$env
    expect_true(.e$impQr)
    .S <- .e$impSamples[[1]]
    .mode1 <- as.numeric(.e$impEtaMode[1, ])
    .H1 <- .e$impEtaHess[[1]]
    .post <- solve(.H1)
    .prop <- .gamma * .post
    # QR accuracy: the sample mean lands on the mode far inside the MC
    # 1/sqrt(N) ~ 0.022 scale
    expect_lt(max(abs(colMeans(.S) - .mode1)), 0.01)
    expect_lt(norm(cov(.S) - .prop, "F") / norm(.prop, "F"), 0.05)
    # the importance weights pull the inflated proposal back to ~H^-1
    .B1 <- .e$impCondVar[[1]]
    expect_equal(unname(.B1), unname(.post), tolerance = 0.02)
    expect_lt(max(abs(.B1 - .post)), max(abs(.B1 - .prop)))
  })

  test_that("Q3: qrShift=FALSE reuses one fixed Sobol point set everywhere", {
    .f <- suppressWarnings(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, "impmap",
              impmapControl(print=0L, nIter=1L, isample=256L, qr=TRUE,
                            qrShift=FALSE)))
    .e <- .f$env
    .z1 <- .zPoints(.e, 1L)
    .z2 <- .zPoints(.e, 2L)
    # identical underlying points across subjects ...
    expect_equal(.z1, .z2, tolerance = 1e-8)
    # ... equal to the raw Sobol N(0,1) point set itself
    expect_equal(unname(.z1), unname(nlmixr2est:::impQrPoints_(256L, 2L, NULL)),
                 tolerance = 1e-8)
  })

  test_that("Q3: qrRefresh pins or redraws the per-subject shift across iterations", {
    .zLast <- function(nIter, qrRefresh, i = 1L) {
      .f <- suppressWarnings(
        nlmixr2(.oneCmt, nlmixr2data::theo_sd, "impmap",
                impmapControl(print=0L, nIter=nIter, isample=128L, qr=TRUE,
                              qrRefresh=qrRefresh)))
      .zPoints(.f$env, i)
    }
    # fixed shift: the same subject's QR points repeat in every iteration
    # (the nIter=2 fit's LAST E-step matches the nIter=1 fit's)
    expect_equal(.zLast(1L, FALSE), .zLast(2L, FALSE), tolerance = 1e-8)
    # refreshed shift: iteration 2 uses different points
    expect_gt(max(abs(.zLast(1L, TRUE) - .zLast(2L, TRUE))), 0.01)
    # even with a fixed shift, different subjects get different shifts
    .f <- suppressWarnings(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, "impmap",
              impmapControl(print=0L, nIter=1L, isample=128L, qr=TRUE,
                            qrRefresh=FALSE)))
    expect_gt(max(abs(.zPoints(.f$env, 1L) - .zPoints(.f$env, 2L))), 0.01)
  })

  test_that("Q3: qr sampler is thread-count independent (D6) and reproducible", {
    .thr0 <- rxode2::getRxThreads()
    on.exit(rxode2::setRxThreads(.thr0), add = TRUE)
    .run <- function(nthr) {
      rxode2::setRxThreads(nthr)
      # perturb the ambient RNG: the fit must pin its own seed (impSeed)
      rxode2::rxSetSeed(sample.int(9999L, 1L)); stats::runif(sample.int(50L, 1L))
      suppressWarnings(
        nlmixr2(.oneCmt, nlmixr2data::theo_sd, "impmap",
                impmapControl(print=0L, nIter=1L, isample=100L,
                              qr=TRUE)))$env$impSamples
    }
    .s1 <- .run(1L)
    .s4 <- .run(4L)
    .s1b <- .run(1L)
    for (.i in seq_along(.s1)) {
      expect_identical(.s1[[.i]], .s4[[.i]])
      expect_identical(.s1[[.i]], .s1b[[.i]])
    }
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
