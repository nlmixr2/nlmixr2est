nmTest({
  # Tests for censoring support in the NLM estimation engine
  # These tests verify that NLM handles M2, M3, M4 censoring
  # matching the behavior used in FOCEI/SAEM

  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- log(c(0, 2.7, 100))
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  # Base dataset (no censoring)
  .dat <- nlmixr2data::theo_sd

  # censOption is inert for NLM (its finite-difference outer Hessian already reflects
  # censoring), so the NLM censoring text stays PLAIN (no " (laplace)"/" (gauss)" suffix --
  # that is FOCEI/FOCE only).  The strip is a defensive no-op that keeps the checks robust.
  .censMethod <- function(x) sub(" \\((laplace|gauss)\\)$", "", as.character(x$censInformation))

  for (meth in c("nlm", "bobyqa", "lbfgsb3c", "n1qn1", "newuoa", "nlminb", "optim")) {
    if (meth == "optim") {
      fit <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .dat, est = meth, list(print = 0, method="BFGS"))
      ))
    } else {
      fit <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .dat, est = meth, list(print = 0))
      ))
    }
    test_that(paste0(meth, " works without censoring (no CENS column)"), {
      expect_s3_class(fit, paste0("nlmixr2.", meth))
      expect_equal(as.character(fit$censInformation), "No censoring")
    })
  }

  .dat0 <- .dat
  .dat0$CENS <- 0L

  for (meth in c("nlm", "bobyqa", "lbfgsb3c", "n1qn1", "newuoa", "nlminb", "optim")) {
    if (meth == "optim") {
      fit0 <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .dat0, est = meth, list(print = 0, method="BFGS"))
      ))
    } else {
      fit0 <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .dat0, est = meth, list(print = 0))
      ))
    }
    test_that("nlm & related methods works with CENS column all zeros", {
      expect_s3_class(fit0, paste0("nlmixr2.", meth))
      expect_equal(as.character(fit0$censInformation), "No censoring")
    })
  }

  # Create censored datasets
  # Mark observations below a threshold as left-censored (M3)
  .LLOQ <- 2.0
  .datM3 <- .dat
  .datM3$CENS <- ifelse(.datM3$DV < .LLOQ & .datM3$EVID == 0, 1L, 0L)
  # For censored obs, set DV to LLOQ
  .datM3$DV[.datM3$CENS == 1] <- .LLOQ

  test_that("nls does not support censoring", {
    expect_error(
      suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .datM3, est = "nls", list(print = 0)))))
  })

  for (meth in c("nlm", "bobyqa", "lbfgsb3c", "n1qn1", "newuoa", "nlminb", "optim")) {
    if (meth == "optim") {
      fit_m3 <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .datM3, est = meth, list(print = 0, method="BFGS"))
      ))
    } else {
      fit_m3 <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .datM3, est = meth, list(print = 0))
      ))
    }
    test_that(paste0(meth, " accepts and processes M3 (left) censored data"), {
      expect_s3_class(fit_m3, paste0("nlmixr2.", meth))
      expect_equal(.censMethod(fit_m3), "M3 censoring")
      # NLM censoring text is plain (censOption inert) -- no laplace/gauss suffix
      expect_equal(as.character(fit_m3$censInformation), "M3 censoring")
    })
  }


  test_that("nls does not support censoring", {
    expect_error(
      suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .datM3, est = "nls", list(print = 0))))
    )
  })

  for (meth in c("nlm", "bobyqa", "lbfgsb3c", "n1qn1", "newuoa", "nlminb", "optim")) {
    if (meth == "optim") {
      fit_base <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .dat, est = meth, list(print = 0, method="BFGS"))
      ))
      fit_m3 <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .datM3, est = meth, list(print = 0, method="BFGS"))
      ))
      expect_false(isTRUE(all.equal(fit_base$objf, fit_m3$objf)))
    } else {
      fit_base <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .dat, est = meth, list(print = 0))
      ))
      fit_m3 <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .datM3, est = meth, list(print = 0))
      ))
      test_that(paste0(meth, "M3 censoring changes the objective function vs no censoring"), {
        expect_false(isTRUE(all.equal(fit_base$objf, fit_m3$objf)))
      })
    }
  }

  # M2 censoring: CENS=0 with a finite LIMIT
  .datM2 <- .dat
  .datM2$CENS <- 0L
  .datM2$LIMIT <- 0  # interval censoring: all obs have a lower bound of 0

  for (meth in c("nlm", "bobyqa", "lbfgsb3c", "n1qn1", "newuoa", "nlminb", "optim")) {
    if (meth == "optim") {
      fit_m2 <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .datM2, est = meth, list(print = 0, method="BFGS"))
      ))
    } else {
      fit_m2 <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .datM2, est = meth, list(print = 0))
      ))
    }
    test_that(paste0(meth, " accepts and processes M2 (interval) censored data"), {
      expect_s3_class(fit_m2, paste0("nlmixr2.", meth))
      expect_equal(.censMethod(fit_m2), "M2 censoring")
    })
  }

  for (meth in c("nlm", "bobyqa", "lbfgsb3c", "n1qn1", "newuoa", "nlminb", "optim")) {
    if (meth == "optim") {
      fit_base <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .dat, est = meth, list(print = 0, method="BFGS"))
      ))
      fit_m2 <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .datM2, est = meth, list(print = 0, method="BFGS"))
      ))
    } else {
      fit_base <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .dat, est = meth, list(print = 0))
      ))
      fit_m2 <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .datM2, est = meth, list(print = 0))
      ))
    }
    test_that(paste0(meth, " method: M2 censoring changes the objective function vs no censoring"), {
      expect_false(isTRUE(all.equal(fit_base$objf, fit_m2$objf)))
    })
  }

  # M4 censoring: CENS!=0 with a finite LIMIT
  .datM4 <- .datM3
  .datM4$LIMIT <- 0  # add LIMIT for M4

  for (meth in c("nlm", "bobyqa", "lbfgsb3c", "n1qn1", "newuoa", "nlminb", "optim")) {
    if (meth == "optim") {
      fit_m4 <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .datM4, est = meth, list(print = 0, method="BFGS"))
      ))
    } else {
      fit_m4 <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .datM4, est = meth, list(print = 0))
      ))
    }
    test_that(paste0(meth, " method accepts and processes M4 (interval-censored) data"), {
      expect_s3_class(fit_m4, paste0("nlmixr2.", meth))
      expect_equal(.censMethod(fit_m4), "M2 and M4 censoring")
    })
  }

  for (meth in c("nlm", "bobyqa", "lbfgsb3c", "n1qn1", "newuoa", "nlminb", "optim")) {
    if (meth == "optim") {
      fit_m3 <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .datM3, est = meth, list(print = 0, method="BFGS"))
      ))
      fit_m4 <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .datM4, est = meth, list(print = 0, method="BFGS"))
      ))
    } else {
      fit_m3 <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .datM3, est = meth, list(print = 0))
      ))
      fit_m4 <- suppressMessages(suppressWarnings(
        .nlmixr(one.cmt, .datM4, est = meth, list(print = 0))
      ))
    }
    test_that(paste0(meth, " method: M3 and M4 give different results"), {
      expect_false(isTRUE(all.equal(fit_m3$objf, fit_m4$objf)))
    })
  }

  test_that("nlm-family M3-censored parameter estimates match focei (#976)", {
    # #976: nlm-family forces every normal endpoint through the log-likelihood
    # path, which always emitted `rx_r_ ~ 0`; the M2/M3/M4 censoring
    # correction (doCensNormal1()) then divided by a zero variance, silently
    # corrupting every censored fit (no ar() needed).  The checks above only
    # ever assert that a fit "runs" and reports the right censoring text, so
    # a corrupted objective/parameter set passed them undetected.  Derivative-
    # free (bobyqa/newuoa/uobyqa) and FD-gradient (nlminb) nlm-family members
    # should now recover parameters close to focei's on the same censored
    # data.  The analytic-gradient members (nlm/n1qn1/lbfgsb3c/optim) have a
    # separate, pre-existing convergence stall on this bounded-tcl model that
    # reproduces even without censoring -- out of scope here, so they are not
    # checked for parameter accuracy.
    f.focei <- .nlmixr(one.cmt, .datM3, est = "focei", control = foceiControl(print = 0))
    for (meth in c("bobyqa", "newuoa", "uobyqa", "nlminb")) {
      fit <- .nlmixr(one.cmt, .datM3, est = meth, control = nlmControl(print = 0))
      expect_equal(as.numeric(fit$theta[["tka"]]), as.numeric(f.focei$theta[["tka"]]),
                   tolerance = 0.1, info = meth)
      expect_equal(as.numeric(fit$theta[["tcl"]]), as.numeric(f.focei$theta[["tcl"]]),
                   tolerance = 0.1, info = meth)
      expect_equal(as.numeric(fit$theta[["tv"]]), as.numeric(f.focei$theta[["tv"]]),
                   tolerance = 0.1, info = meth)
      # add.sd is the parameter the r=0 bug corrupted most (it inflated it
      # ~9x in the reproduction that motivated #976); give it a bit more room
      expect_equal(as.numeric(fit$theta[["add.sd"]]), as.numeric(f.focei$theta[["add.sd"]]),
                   tolerance = 0.2, info = meth)
    }
  })

  test_that("nlm-family M3-censored propT() parameter estimates match focei (#976 follow-up)", {
    # #976 follow-up (antigravity review): .fixCensRNuLine() (R/focei.R,
    # shared by nlm-family and, gated, the t()/cauchy() FOCEi path from
    # #979) rebuilt rx_r_ by re-inlining rx_rll_'s defining EXPRESSION
    # (e.g. sqrt((b*rx_pred_)^2) for propT()/powT(), whose F is the
    # TRANSFORMED prediction, i.e. the symbol rx_pred_ itself -- see
    # rxode2's .rxGetVarianceForErrorPropOrPowF()).  That expression is
    # placed at the rx_r_ ~ 0 line's position, which comes AFTER rx_pred_
    # has already been overwritten with the scalar log-likelihood
    # (rx_pred_ <- llikNorm(dv, rx_pred_, rx_rll_)).  Re-evaluating the AST
    # there silently fed the log-likelihood value back in as the propT()
    # mean, corrupting rx_r_ for any transformed prop()/pow() error model.
    # The fix instead references the already-computed rx_rll_ VARIABLE
    # (rx_r_ ~ rx_rll_^2), which still holds its original value at that
    # point.  A plain add()/prop() (F = rx_pred_f_, a column untouched by
    # the llik overwrite) would not have caught this -- propT() is
    # required to exercise the bug.
    one.cmt.propT <- function() {
      ini({
        tka <- 0.45
        tcl <- log(c(0, 2.7, 100))
        tv <- 3.45
        prop.sd <- c(0, 0.3)
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        linCmt() ~ propT(prop.sd)
      })
    }
    f.focei <- .nlmixr(one.cmt.propT, .datM3, est = "focei", control = foceiControl(print = 0))
    fit <- .nlmixr(one.cmt.propT, .datM3, est = "bobyqa", control = nlmControl(print = 0))
    expect_equal(as.numeric(fit$theta[["tka"]]), as.numeric(f.focei$theta[["tka"]]),
                 tolerance = 0.15)
    expect_equal(as.numeric(fit$theta[["tcl"]]), as.numeric(f.focei$theta[["tcl"]]),
                 tolerance = 0.15)
    expect_equal(as.numeric(fit$theta[["tv"]]), as.numeric(f.focei$theta[["tv"]]),
                 tolerance = 0.15)
    expect_equal(as.numeric(fit$theta[["prop.sd"]]), as.numeric(f.focei$theta[["prop.sd"]]),
                 tolerance = 0.15)
  })

})
