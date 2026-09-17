# est="ifocei" and est="mfocei" were the only two FOCEi family methods without
# the "iov" attribute, so .uiApplyIov() stood down and nothing expanded the
# occasion parameters -- every IOV model errored under them (#1083).
test_that("IOV models fit with ifocei and mfocei (#1083)", {
  skip_on_cran()

  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.3
      iov.ka ~ 0.1 | occ
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka + iov.ka)
      cl <- exp(tcl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  theoIov <- nlmixr2data::theo_sd
  theoIov$occ <- 1L + (theoIov$TIME >= 5)

  for (.est in c("ifocei", "mfocei")) {
    .fit <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, theoIov, est = .est,
              control = list(print = 0L, maxOuterIterations = 0L,
                             covMethod = "", calcTables = FALSE))))
    # calcTables=FALSE, so the fit is the core object rather than the data frame
    expect_s3_class(.fit, "nlmixr2FitCore")
    # the occasion parameter is expanded for the fit and restored afterwards
    expect_true("iov.ka" %in% .fit$ui$iniDf$name, info = .est)
    expect_equal(.fit$ui$iniDf$condition[.fit$ui$iniDf$name == "iov.ka"], "occ",
                 info = .est)
  }
})
