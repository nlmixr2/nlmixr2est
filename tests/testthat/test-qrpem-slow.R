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

  test_that("impSeed decorrelates the draws yet stays reproducible", {
    .draw <- function(seed) {
      suppressWarnings(
        nlmixr2(.oneCmt, nlmixr2data::theo_sd, "impmap",
                impmapControl(print=0L, nIter=1L, isample=100L,
                              impSeed=seed)))$env$impSamples[[1]]
    }
    # a fixed impSeed is bit-reproducible ...
    expect_identical(.draw(1L), .draw(1L))
    # ... and a different impSeed gives a genuinely different sample set
    # (the pre-fix behavior ignored impSeed entirely -- samples were identical)
    expect_false(identical(.draw(1L), .draw(999L)))
    expect_gt(max(abs(.draw(1L) - .draw(999L))), 0.05)
  })

  test_that("Q4: impCov=TRUE with qr=TRUE gives an SPD covariance matching FOCEI |r|", {
    .fi <- suppressWarnings(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, "impmap",
              impmapControl(print=0L, nIter=40L, isample=500L, qr=TRUE,
                            impCov=TRUE)))
    .se <- as.numeric(.fi$env$impSe)
    .nth <- .fi$env$impCovThetaN
    expect_true(all(is.finite(.se) & .se > 0))
    .cov <- .fi$env$impCov
    expect_equal(.cov, t(.cov), tolerance = 1e-8)
    expect_true(all(eigen(.cov, symmetric = TRUE, only.values = TRUE)$values > 0))
    .ff <- suppressWarnings(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, "focei",
              foceiControl(print = 0L, covMethod = "r")))
    skip_if(is.null(.ff$cov), "FOCEI |r| covariance unavailable")
    .fse <- sqrt(diag(.ff$cov))[seq_len(.nth)]
    expect_equal(.se[seq_len(.nth)], unname(.fse), tolerance = 0.1)
  })

  test_that("Q6: sirSample = isample reproduces the full weighted M-step exactly", {
    # the SIR branch only engages when sirN < isample; at equality the fit is
    # bit-identical to sir=FALSE (the SIR plumbing is a strict superset)
    .run <- function(sir, sirSample = NULL) {
      .ctl <- if (sir) {
        impmapControl(print=0L, nIter=2L, isample=100L, sir=TRUE,
                      sirSample=sirSample)
      } else {
        impmapControl(print=0L, nIter=2L, isample=100L)
      }
      suppressWarnings(nlmixr2(.oneCmt, nlmixr2data::theo_sd, "impmap", .ctl))
    }
    .f0 <- .run(FALSE)
    .fEq <- .run(TRUE, 100L)
    expect_equal(fixef(.f0), fixef(.fEq), tolerance = 1e-10)
    expect_equal(.f0$env$impObj, .fEq$env$impObj, tolerance = 1e-10)
    # a genuine resample (half the samples) stays close after one M-step
    .fH <- .run(TRUE, 50L)
    expect_lt(max(abs(fixef(.fH) - fixef(.f0))), 0.05)
  })

  test_that("Q6: sir=TRUE converges to FOCEI (non-mu theta + sigma unfixed)", {
    .mstr <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp ~ add(add.sd)
      })
    }
    .d <- nlmixr2data::theo_sd
    .ff <- suppressWarnings(nlmixr2(.mstr, .d, "focei",
                                    foceiControl(print = 0L, covMethod = "")))
    # SIR with the default resample size (30 of 300)
    .fs <- suppressWarnings(nlmixr2(.mstr, .d, "impmap",
                                    impmapControl(print = 0L, nIter = 30L,
                                                  isample = 300L, sir = TRUE)))
    expect_true(.fs$env$impSir)
    expect_identical(.fs$env$impSirSample, 30L)
    expect_equal(unname(fixef(.fs)["tv"]), unname(fixef(.ff)["tv"]), tolerance = 0.05)
    expect_equal(unname(fixef(.fs)["add.sd"]), unname(fixef(.ff)["add.sd"]), tolerance = 0.05)
    expect_equal(fixef(.fs)[c("tka", "tcl")], fixef(.ff)[c("tka", "tcl")], tolerance = 0.05)
    expect_equal(unname(diag(.fs$omega)), unname(diag(.ff$omega)), tolerance = 0.1)
    # qr + sir combined (the full QRPEM configuration)
    .fq <- suppressWarnings(nlmixr2(.mstr, .d, "impmap",
                                    impmapControl(print = 0L, nIter = 30L,
                                                  isample = 300L, qr = TRUE,
                                                  sir = TRUE)))
    expect_equal(unname(fixef(.fq)["tv"]), unname(fixef(.ff)["tv"]), tolerance = 0.05)
    expect_equal(unname(fixef(.fq)["add.sd"]), unname(fixef(.ff)["add.sd"]), tolerance = 0.05)
    expect_equal(fixef(.fq)[c("tka", "tcl")], fixef(.ff)[c("tka", "tcl")], tolerance = 0.05)
  })

  test_that("Q7: est='qrpem' equals impmap(qr+sir) and matches FOCEI", {
    .fq <- suppressWarnings(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, "qrpem",
              qrpemControl(print = 0L, nIter = 30L, isample = 300L)))
    expect_true(.fq$env$impQr)
    expect_true(.fq$env$impSir)
    expect_identical(.fq$env$method, "qrpem")
    # the sugar is exactly impmap with qr=TRUE, sir=TRUE (same seed pinning)
    .fi <- suppressWarnings(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, "impmap",
              impmapControl(print = 0L, nIter = 30L, isample = 300L,
                            qr = TRUE, sir = TRUE)))
    expect_equal(fixef(.fq), fixef(.fi), tolerance = 1e-10)
    expect_equal(.fq$env$impObj, .fi$env$impObj, tolerance = 1e-10)
    # and it lands on the FOCEI estimates
    .ff <- suppressWarnings(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, "focei",
              foceiControl(print = 0L, covMethod = "")))
    expect_equal(fixef(.fq), fixef(.ff), tolerance = 0.05)
    expect_equal(unname(diag(.fq$omega)), unname(diag(.ff$omega)), tolerance = 0.1)
  })

  test_that("Q1: default (qr/sir off) fit is unchanged vs the pre-QRPEM baseline", {
    .ref <- readRDS(test_path("baselines", "qrpem-baseline-ref.rds"))
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
