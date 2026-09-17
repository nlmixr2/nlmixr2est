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

  .fitIt <- function(est) {
    suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, theoIov, est = est,
              control = list(print = 0L, maxOuterIterations = 0L,
                             covMethod = "", calcTables = FALSE))))
  }

  # the "*f" delegate is the SAME method with foceiControl(fast=TRUE): it only
  # changes how the outer gradient is computed, so at a fixed starting point the
  # objective must agree exactly.  That is the assertion that the expansion is
  # right -- the restored iniDf shape below would look fine either way, and
  # ifoceif/mfoceif already fitted this model before the fix.
  for (.est in c("ifocei", "mfocei")) {
    .fit <- .fitIt(.est)
    # calcTables=FALSE, so the fit is the core object rather than the data frame
    expect_s3_class(.fit, "nlmixr2FitCore")
    expect_equal(.fit$objf, .fitIt(paste0(.est, "f"))$objf, info = .est)
    # the occasion parameter is expanded for the fit and restored afterwards
    expect_true("iov.ka" %in% .fit$ui$iniDf$name, info = .est)
    expect_equal(.fit$ui$iniDf$condition[.fit$ui$iniDf$name == "iov.ka"], "occ",
                 info = .est)
    # and the model line is the user's again, with no rewrite residue
    .txt <- paste(vapply(.fit$ui$lstExpr,
                         function(x) paste(deparse(x), collapse = " "),
                         character(1)), collapse = " ")
    expect_false(grepl("rx.iov.", .txt, fixed = TRUE), info = .est)
  }
})

# The same registry drift one file over: the "full" conditional-Hessian
# delegates refused a correlated occasion block that the base method they
# dispatch to fits (#1083).
test_that("a correlated occasion block fits under a full Laplace/AGQ delegate", {
  skip_on_cran()

  corr.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      iov.cl + iov.v ~ c(0.1, 0.03, 0.2) | occ
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + iov.cl)
      v <- exp(tv + iov.v)
      linCmt() ~ add(add.sd)
    })
  }

  theoIov <- nlmixr2data::theo_sd
  theoIov$occ <- 1L + (theoIov$TIME >= 5)

  for (.est in c("laplace", "flaplace", "fagq")) {
    .fit <- suppressMessages(suppressWarnings(
      nlmixr2(corr.cmt, theoIov, est = .est,
              control = list(print = 0L, maxOuterIterations = 0L,
                             covMethod = "", calcTables = FALSE))))
    expect_s3_class(.fit, "nlmixr2FitCore")
    expect_true(is.finite(.fit$objf), info = .est)
    # both occasion parameters come back on their own `| occ` rows
    expect_equal(.fit$ui$iniDf$condition[.fit$ui$iniDf$name %in%
                                           c("iov.cl", "iov.v")],
                 c("occ", "occ"), info = .est)
  }
})
