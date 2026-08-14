# foceif must reach the SAME optimum as focei, whichever outer optimizer it uses.
#
# It did not: lbfgsb3c ran with the projected-gradient test suppressed, leaving
# `lbfgsFactr` -- a single-step objective-reduction test, not a stationarity test
# -- as its only stopping rule, so the fit stopped after roughly as many outer
# evaluations as the model has parameters and reported a worse objective.  That
# is enough to reorder a set of candidate models (discussion #924).
#
# Fit-based, so weekly-batched; the control-level checks are in
# test-focei-pgtol-control.R (push/PR subset).

nmTest({
  .oneCmt <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  # counts one row per outer objective evaluation
  .outerEvals <- function(fit) sum(fit$parHistData$type == "Unscaled")
  .gradTypes <- function(fit) {
    unique(as.character(fit$parHistData$type[
      grepl("Gradient|Difference|Sensitivity", fit$parHistData$type)]))
  }

  test_that("foceif converges to focei's optimum under either outer optimizer", {
    .focei <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(print = 0))))
    .lbfgs <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, est = "foceif",
              control = foceiControl(print = 0))))
    .nlminb <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, est = "foceif",
              control = foceiControl(print = 0, outerOpt = "nlminb"))))

    # the mechanism actually ran -- the counter/type proves the analytic gradient
    # was used, which the OFV agreement alone would not
    expect_true(all(grepl("^Analytic Gradient", .gradTypes(.lbfgs))))
    expect_true(all(grepl("^Analytic Gradient", .gradTypes(.nlminb))))
    # ... and that the outer projected-gradient test was actually armed
    expect_true(.lbfgs$control$lbfgsPgtolOuter > 0)

    # no outer optimizer may be more than ~1 OFV unit worse than the best of them;
    # the defect was worth 0.3-1.1 units, which reorders candidate models
    .best <- min(.focei$objf, .lbfgs$objf, .nlminb$objf)
    expect_lt(.focei$objf - .best, 1)
    expect_lt(.lbfgs$objf - .best, 1)
    expect_lt(.nlminb$objf - .best, 1)

    # the old failure signature: stopping at about one evaluation per parameter.
    # 7 parameters are estimated here (4 thetas + 3 omegas).
    expect_gt(.outerEvals(.lbfgs), 7)
  })
})
