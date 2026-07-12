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
