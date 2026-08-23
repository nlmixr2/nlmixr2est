# linCmt() sensitivity carry with event-modifier (f()/alag()) channels and a
# non-cl/v parameterization, validated at the fit level against rxode2's own
# linToOde() ODE translation integrated with useLinCmt=FALSE (the ODE side
# carries f()/alag() as genuine event sensitivities and every trans exactly),
# plus the data-aware fallbacks the jump channels need.  Gradient-level checks
# are in test-focei-lincmt-carry-jump.R; shared fixtures in
# helper-lincmt-carry.R.
#
# Slow batch (real FOCEi fits) -- see .slowBatches in tests/testthat.R.
# Everything here needs an rxode2 with the carry sentinels and the jump pin
# and skips cleanly on a released rxode2 without them.

.carryJumpFitCapable <- function() {
  .rxFoceiLinCmtCarryCapable() && .rxFoceiLinCmtCarryJumpCapable() # nolint: object_usage_linter.
}

.carryJumpFit <- function(ui, dat, carry, maxOut = 0L) {
  suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(ui, dat, est = "focei",
                        control = .carryFitCtl(carry, maxOut)))) # nolint: object_usage_linter.
}

# simulate DV from an explicitly integrated ODE truth under nocb
.carryJumpSimDv <- function(odeTxt, params, dat, nid, sd = 0.3) {
  m <- rxode2::rxode2(odeTxt)
  dv <- unlist(lapply(seq_len(nid), function(i) {
    rxode2::rxSolve(m, params = params(i), events = dat[dat$id == i, ],
                    returnType = "data.frame", covsInterpolation = "nocb",
                    useLinCmt = FALSE, atol = 1e-10, rtol = 1e-10)$cp
  }))
  obs <- dat$evid == 0
  dat$dv <- 0
  set.seed(99)
  dat$dv[obs] <- dv + rnorm(sum(obs), 0, sd)
  dat
}

test_that("f()+alag()+covariate fit matches the ODE reference only with the carry", {
  skip_on_cran()
  skip_if_not(.carryJumpFitCapable())
  dat <- .carryFitDat()
  dat$cmt <- 1
  set.seed(17)
  et <- matrix(rnorm(18, 0, 0.3), 6)
  dat <- .carryJumpSimDv("
cl = exp(tcl)*(wt/70)^0.75*exp(eta_cl)
v = exp(tv)
ka = exp(tka)
f(depot) = expit(tf + eta_f + (wt-70)/70)
alag(depot) = exp(tlag + eta_lag)
d/dt(depot) = -ka*depot
d/dt(central) = ka*depot - (cl/v)*central
cp = central/v", function(i) {
    c(tcl = log(2), tv = log(20), tka = log(1.2), tf = -0.5, tlag = log(0.5),
      eta_cl = et[i, 1], eta_f = et[i, 2], eta_lag = et[i, 3])
  }, dat, 6L)
  uiO <- rxode2::linToOde(rxode2::rxode2(.carryModJump))
  fO <- .carryJumpFit(uiO, dat, "none")
  fC <- .carryJumpFit(.carryModJump, dat, "auto")
  fN <- .carryJumpFit(.carryModJump, dat, "none")
  expect_true(grepl("rx_lcCarryPin_", paste(rxode2::rxNorm(fC$env$innerModel), collapse = "")))
  # measured: carry 4.9e-5 vs naive 0.127
  expect_lt(abs(fC$objective - fO$objective), 0.01)
  expect_gt(abs(fN$objective - fO$objective), 0.05)
  etaC <- as.matrix(fC$eta[, -1])
  etaO <- as.matrix(fO$eta[, -1])
  expect_lt(max(abs(etaC - etaO)), 2e-3)
  # full fits: this six-subject lag+F problem is weakly identified and the
  # inner problem is multi-modal in eta.lag, so the two outer runs end in
  # different basins (warm-started etas along different paths).  The
  # statement that holds is surface agreement: a fresh posthoc evaluation
  # of both models at EITHER optimum gives the same objective (measured
  # 3e-4 apart; naive 0.14-0.17 apart).
  FO <- .carryJumpFit(uiO, dat, "none", 200L)
  FC <- .carryJumpFit(.carryModJump, dat, "auto", 200L)
  setTh <- function(ui, th) do.call(rxode2::ini, c(list(ui), as.list(th)))
  uiC <- rxode2::rxode2(.carryModJump)
  for (th in list(FO$theta, FC$theta)) {
    pO <- .carryJumpFit(setTh(uiO, th), dat, "none")
    pC <- .carryJumpFit(setTh(uiC, th), dat, "auto")
    pN <- .carryJumpFit(setTh(uiC, th), dat, "none")
    expect_lt(abs(pC$objective - pO$objective), 0.01)
    expect_gt(abs(pN$objective - pO$objective), 0.05)
  }
})

test_that("a 2-cmt A/B/alpha/beta model with a covariate on B fits like its ODE", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  mod <- function() {
    ini({ta <- log(0.5); tb <- log(0.02); tA <- log(0.04); tB <- log(0.01)
         eta.B ~ 0.1; add.sd <- 0.3})
    model({
      alpha <- exp(ta); beta <- exp(tb); A <- exp(tA)
      B <- exp(tB) * (wt / 70)^-1 * exp(eta.B)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  dat <- .carryFitDat()
  set.seed(23)
  et <- rnorm(6, 0, 0.3)
  uiO <- rxode2::linToOde(rxode2::rxode2(mod))
  odeTxt <- paste(rxode2::rxNorm(uiO), collapse = "\n")
  odeTxt <- gsub("eta.B", "eta_B", odeTxt, fixed = TRUE)
  dat <- .carryJumpSimDv(odeTxt, function(i) {
    c(ta = log(0.5), tb = log(0.02), tA = log(0.04), tB = log(0.01), eta_B = et[i])
  }, dat, 6L, sd = 0.1)
  fO <- .carryJumpFit(uiO, dat, "none")
  fC <- .carryJumpFit(mod, dat, "auto")
  fN <- .carryJumpFit(mod, dat, "none")
  expect_true(grepl("rx_lcCarryAdv_", paste(rxode2::rxNorm(fC$env$innerModel), collapse = "")))
  expect_lt(abs(fC$objective - fO$objective), 0.01)
  expect_gt(abs(fN$objective - fO$objective), abs(fC$objective - fO$objective) * 10)
})

test_that("doses outside the f()/alag() compartment fall back with a runInfo note", {
  skip_on_cran()
  skip_if_not(.carryJumpFitCapable())
  d <- .carryFitDat(2L)
  d$cmt <- ifelse(d$evid == 1 & d$time > 20, 2, 1) # later doses go to central
  d$dv <- ifelse(d$evid == 0, 3, 0)
  fit <- .carryJumpFit(.carryModJump, d, "auto")
  expect_identical(fit$foceiControl$linCmtSensCarry, "none")
  expect_true(any(grepl("carry gradient off", unlist(fit$runInfo))))
})

test_that("an infusion with a lag channel falls back with a runInfo note", {
  skip_on_cran()
  skip_if_not(.carryJumpFitCapable())
  d <- .carryFitDat(2L)
  d$cmt <- 1
  d$rate <- ifelse(d$evid == 1, d$amt / 2, 0)
  d$dv <- ifelse(d$evid == 0, 3, 0)
  fit <- .carryJumpFit(.carryModJump, d, "auto")
  expect_identical(fit$foceiControl$linCmtSensCarry, "none")
  expect_true(any(grepl("carry gradient off", unlist(fit$runInfo))))
})
