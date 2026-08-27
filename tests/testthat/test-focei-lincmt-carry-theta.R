test_that("theta-side carry: eligibility, emission and FD (#1003)", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())

  modTheta <- function() {
    ini({
      tcl <- log(2)
      tv <- log(20)
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl) * (wt / 70)^0.75 # nolint: object_usage_linter.
      v <- exp(tv) # nolint: object_usage_linter.
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  modThetaNoCov <- function() {
    ini({
      tcl <- log(2)
      tv <- log(20)
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl) # nolint: object_usage_linter.
      v <- exp(tv) # nolint: object_usage_linter.
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }

  buildTxt <- function(mod, carry, method = "nlm") {
    ui <- suppressMessages(nlmixr2est::nlmixr2(mod))
    u <- rxode2::.copyUi(ui)
    ctl <- if (method == "nlm") {
      nlmixr2est::nlmControl(linCmtSensCarry = carry)
    } else {
      nlmixr2est::nlsControl(linCmtSensCarry = carry)
    }
    assign("control", ctl, envir = u)
    e <- suppressMessages(if (method == "nlm") u$nlmEnv else u$nlsEnv)
    if (method == "nlm") e$..nlmS else e$..nlsS
  }

  # mechanism: the covariate theta's line is substituted, and only then
  txtC <- buildTxt(modTheta, "auto")
  txtN <- buildTxt(modTheta, "none")
  expect_true(any(grepl("rx_lcConc_~linCmtB", strsplit(txtC, "\n")[[1]])))
  expect_true(any(grepl("rx_lcCarryS0r0_", strsplit(txtC, "\n")[[1]])))
  expect_false(any(grepl("rx_lcCarry|rx_lcConc_", strsplit(txtN, "\n")[[1]])))

  # a model with no covariate on any slot is byte-identical with carry on/off
  expect_identical(buildTxt(modThetaNoCov, "auto"), buildTxt(modThetaNoCov, "none"))

  # FD on theta through the generated gradient model
  ev <- .carryEv()
  ev$dv <- ifelse(ev$evid == 0, c(0, 2, 0, 1.5, 0, 1.2, 0, 0.9)[seq_len(nrow(ev))], 0)
  pars <- c(`THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = 0.5)
  fdErr <- function(txt) {
    m <- suppressWarnings(rxode2::rxode2(txt))
    slv <- function(q) {
      rxode2::rxSolve(m,
        params = q, events = ev, returnType = "data.frame",
        covsInterpolation = "nocb"
      )
    }
    r0 <- slv(pars)
    h <- 1e-5
    vapply(1:3, function(k) {
      got <- r0[[paste0("rx__sens_rx_pred__BY_THETA_", k, "___")]]
      if (is.null(got)) {
        return(0)
      } # nls carries no residual-sd column
      tn <- paste0("THETA[", k, "]")
      a <- pars
      a[tn] <- a[tn] + h
      b <- pars
      b[tn] <- b[tn] - h
      fd <- (slv(a)$rx_pred_ - slv(b)$rx_pred_) / (2 * h)
      max(abs(got - fd) / (abs(fd) + 1e-8))
    }, numeric(1))
  }
  errC <- fdErr(txtC)
  errN <- fdErr(txtN)
  expect_true(all(errC < 1e-6))
  # the naive covariate-theta gradient is visibly wrong on the same data
  expect_true(errN[1] > 1e-3)

  # nls: same substitution and FD for the covariate theta
  txtCnls <- buildTxt(modTheta, "auto", "nls")
  expect_true(any(grepl("rx_lcConc_~linCmtB", strsplit(txtCnls, "\n")[[1]])))
  errCnls <- fdErr(txtCnls)
  expect_true(all(errCnls[1:2] < 1e-6))
})

test_that("theta-side carry eligibility keeps bias-to-false rules", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())

  # a non-separable theta (estimated allometric exponent) is rejected
  modExp <- function() {
    ini({
      tcl <- log(2)
      texp <- 0.75
      tv <- log(20)
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl) * (wt / 70)^texp # nolint: object_usage_linter.
      v <- exp(tv) # nolint: object_usage_linter.
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  ui <- suppressMessages(nlmixr2est::nlmixr2(modExp))
  s <- suppressMessages(ui$nlmThetaS)
  thetaVars <- paste0("THETA_", seq_len(s$..maxTheta), "_")
  fp <- .rxCarryFactorPred(s)
  expect_false(is.null(fp))
  pairs <- .rxCarryThetaEligible(list(ui), s, thetaVars, fp)
  # two estimated thetas share the cl slot, so BOTH are rejected (one
  # direction per slot; the non-separable texp could not qualify anyway)
  expect_identical(nrow(pairs), 0L)

  # fixing the exponent leaves tcl alone in the slot: eligible again
  modExpFix <- function() {
    ini({
      tcl <- log(2)
      texp <- fix(0.75)
      tv <- log(20)
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl) * (wt / 70)^texp # nolint: object_usage_linter.
      v <- exp(tv) # nolint: object_usage_linter.
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  ui2 <- suppressMessages(nlmixr2est::nlmixr2(modExpFix))
  s2 <- suppressMessages(ui2$nlmThetaS)
  thetaVars2 <- paste0("THETA_", seq_len(s2$..maxTheta), "_")
  fp2 <- .rxCarryFactorPred(s2)
  pairs2 <- .rxCarryThetaEligible(list(ui2), s2, thetaVars2, fp2)
  expect_true("tcl" %in% pairs2$etaName)
})
