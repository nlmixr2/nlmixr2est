# issue 1004, fit level: an ll() linCmt() endpoint under a time-varying covariate
# against the same ll() endpoint on the linToOde() ODE equivalent,
# integrated with useLinCmt=FALSE (genuine sensitivities).  Slow batch -- see
# .slowBatches in tests/testthat.R; shared fixtures in helper-lincmt-carry.R.

test_that("ll() carry fit matches the ODE reference; naive does not", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  dat <- .carryFitDat()
  uiL <- suppressMessages(nlmixr2est::nlmixr2(.carryModLlNorm))
  uiO <- rxode2::linToOde(uiL)
  # the translation must keep the explicit likelihood line
  expect_true(any(grepl("^ll\\(", as.character(uiO$lstExpr))))
  m <- rxode2::rxode2("
cl = exp(tcl)*(wt/70)^0.75*exp(eta_cl)
v = exp(tv)
d/dt(central) = -(cl/v)*central
cp = central/v")
  set.seed(17)
  etaTrue <- rnorm(6, 0, 0.3)
  dv <- unlist(lapply(1:6, function(i) {
    rxode2::rxSolve(m, params = c(tcl = log(2), tv = log(20),
                                  eta_cl = etaTrue[i]),
                    events = dat[dat$id == i, ], returnType = "data.frame",
                    covsInterpolation = "nocb", useLinCmt = FALSE)$cp
  }))
  obs <- dat$evid == 0
  dat$dv <- 0
  set.seed(99)
  dat$dv[obs] <- dv + rnorm(sum(obs), 0, 0.3)
  fit <- function(ui, carry, maxOut = 0L) {
    suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(ui, dat, est = "focei",
                          control = .carryFitCtl(carry, maxOut))))
  }
  fO <- fit(uiO, "none")
  fC <- fit(uiL, "auto")
  fN <- fit(uiL, "none")
  expect_true(grepl("rx_lcConc_",
                    paste(rxode2::rxNorm(fC$env$innerModel), collapse = "\n")))
  gapC <- abs(fC$objective - fO$objective)
  gapN <- abs(fN$objective - fO$objective)
  expect_lt(gapC, 0.01)
  expect_gt(gapN, 0.05)
  expect_lt(max(abs(fC$eta$eta.cl - fO$eta$eta.cl)), 1e-3)
  fo <- fit(uiO, "none", 200L)
  fc <- fit(uiL, "auto", 200L)
  expect_lt(abs(FC$objective - FO$objective), 0.01)
  expect_lt(max(abs(FC$theta - FO$theta)), 5e-3)
})
