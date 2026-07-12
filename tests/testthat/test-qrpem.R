# Quick QRPEM core tests (no model fits) -- always run on push/PR CI.
# Fit-based QRPEM tests live in test-qrpem-slow.R (weekly slow-tests batch).
nmTest({
  test_that("impmapControl qr/sir option defaults and validation", {
    .ctl <- impmapControl()
    expect_false(.ctl$qr)
    expect_true(.ctl$qrShift)
    expect_true(.ctl$qrRefresh)
    expect_false(.ctl$sir)
    # sirSample NULL resolves to max(25, ceiling(isample/10)) at build time
    expect_identical(.ctl$sirSample, 30L)
    expect_identical(impmapControl(isample=100L)$sirSample, 25L)
    expect_identical(impmapControl(isample=4000L)$sirSample, 400L)
    expect_identical(impmapControl(sirSample=50L)$sirSample, 50L)

    expect_error(impmapControl(qr="yes"))
    expect_error(impmapControl(qrShift=NA))
    expect_error(impmapControl(qrRefresh=c(TRUE, FALSE)))
    expect_error(impmapControl(sir=1))
    expect_error(impmapControl(sirSample=0L))
    expect_error(impmapControl(sirSample=301L), "cannot exceed")

    # round-trips through do.call (getValidNlmixrCtl / .foceiFamilyControl)
    .ctl <- impmapControl(qr=TRUE, qrShift=FALSE, qrRefresh=FALSE,
                          sir=TRUE, sirSample=40L)
    .ctl2 <- do.call(impmapControl, .ctl)
    expect_true(.ctl2$qr)
    expect_false(.ctl2$qrShift)
    expect_false(.ctl2$qrRefresh)
    expect_true(.ctl2$sir)
    expect_identical(.ctl2$sirSample, 40L)

    # impControl() inherits the options through ...
    .ic <- impControl(qr=TRUE, sir=TRUE)
    expect_true(.ic$qr)
    expect_true(.ic$sir)
    expect_identical(.ic$mapIter, 0L)
  })

  test_that("impQrPoints_ produces low-discrepancy N(0,1) points", {
    .Z <- nlmixr2est:::impQrPoints_(256L, 2L, NULL)
    expect_true(is.matrix(.Z) && all(dim(.Z) == c(256L, 2L)))
    expect_true(all(is.finite(.Z)))
    # deterministic: same call, same points
    expect_identical(.Z, nlmixr2est:::impQrPoints_(256L, 2L, NULL))
    # low-discrepancy signature: qnorm(sobol) column means are O(log N / N),
    # far below the 1/sqrt(N) = 0.0625 pseudo-random scale
    expect_true(max(abs(colMeans(.Z))) < 0.01)
    # near-perfect stratification: the boost engine skips the zero point, so a
    # 2^8 block is offset by one -- 16 equal bins hold 16 +/- 2 points each
    # (a 256-point pseudo-random draw would routinely miss by 8+)
    .U <- pnorm(.Z)
    for (.j in 1:2) {
      .cnt <- table(cut(.U[, .j], breaks = seq(0, 1, by = 1/16)))
      expect_true(all(abs(.cnt - 16L) <= 2L))
    }
    # unit variance to QR accuracy
    expect_equal(unname(apply(.Z, 2, sd)), c(1, 1), tolerance = 0.02)
  })

  test_that("impQrPoints_ Cranley-Patterson shift wraps and stays stratified", {
    .Z0 <- nlmixr2est:::impQrPoints_(256L, 2L, NULL)
    .sh <- c(0.371, 0.842)
    .Z <- nlmixr2est:::impQrPoints_(256L, 2L, .sh)
    expect_true(all(is.finite(.Z)))
    expect_false(identical(.Z, .Z0))
    # the shift acts mod 1 on the uniforms
    .U0 <- pnorm(.Z0)
    .U <- pnorm(.Z)
    expect_equal(.U, (.U0 + rep(.sh, each = 256L)) %% 1, tolerance = 1e-8)
    # shifted Sobol keeps near-perfect stratification (16 +/- 2 per 16 bins)
    for (.j in 1:2) {
      .cnt <- table(cut(.U[, .j], breaks = seq(0, 1, by = 1/16)))
      expect_true(all(abs(.cnt - 16L) <= 2L))
    }
    # a shift near 1 wraps rather than escaping (0,1)
    .Zw <- nlmixr2est:::impQrPoints_(64L, 2L, c(0.999999, 0.5))
    expect_true(all(is.finite(.Zw)))
    # bad input
    expect_error(nlmixr2est:::impQrPoints_(256L, 2L, c(0.5)))
    expect_error(nlmixr2est:::impQrPoints_(0L, 2L, NULL))
  })

  test_that("qr/sir names are stripped when down-converting to foceiControl", {
    .env <- new.env()
    .env$impmapControl <- impmapControl(qr=TRUE, sir=TRUE)
    .fc <- nlmixr2est:::.impmapControlToFoceiControl(.env, assign=FALSE)
    expect_s3_class(.fc, "foceiControl")
    for (.n in c("qr", "qrShift", "qrRefresh", "sir", "sirSample")) {
      expect_null(.fc[[.n]])
    }
  })
})
