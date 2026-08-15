# foceif must reach focei's optimum, not stop 8.6 OFV units short of it.
#
# A 2-cmt oral fit is the case that pins this: at the old lbfgsFactr default the
# analytic-gradient path stopped after 6 outer evaluations, +8.63 OFV above the
# bobyqa reference, which is enough to reorder candidate models (discussion
# #924).  Two orders tighter it reaches -0.10 in 13 evaluations.
#
# Fit-based, so weekly-batched; the control-level check is in
# test-focei-factr-control.R (push/PR subset).

nmTest({
  .twoCmtOral <- function() {
    ini({
      tka <- 0.5; tcl <- 1.0; tv <- 4.0; tvp <- 4.0; tq <- 1.0
      eta.ka ~ 0.3; eta.cl ~ 0.3; eta.v ~ 0.1
      prop.sd <- 0.2
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      vp <- exp(tvp); q <- exp(tq)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center - q / v * center + q / vp * peri
      d/dt(peri)   <-  q / v * center - q / vp * peri
      cp <- center / v
      cp ~ prop(prop.sd)
    })
  }
  .dat <- nlmixr2data::Oral_2CPT[nlmixr2data::Oral_2CPT$SD == 1, ]

  .outerEvals <- function(fit) sum(fit$parHistData$type == "Unscaled")
  .gradTypes <- function(fit) {
    unique(as.character(fit$parHistData$type[
      grepl("Gradient|Difference|Sensitivity", fit$parHistData$type)]))
  }

  test_that("foceif reaches focei's optimum on a 2-cmt oral fit", {
    # sigdig is pinned because the absolute objective below is a numeric literal
    # and sigdig drives both the optimizer tolerances and the ODE rtol/atol
    .focei <- suppressWarnings(suppressMessages(
      nlmixr2(.twoCmtOral, .dat, est = "focei",
              control = foceiControl(print = 0, sigdig = 3))))
    .foceif <- suppressWarnings(suppressMessages(
      nlmixr2(.twoCmtOral, .dat, est = "foceif",
              control = foceiControl(print = 0, sigdig = 3))))

    # the mechanism ran -- OFV agreement alone would not show the analytic
    # gradient was ever used
    expect_true(all(grepl("^Analytic Gradient", .gradTypes(.foceif))))

    # the defect was +8.63; anything under 1 OFV unit is not going to reorder
    # candidate models
    expect_lt(.foceif$objf - .focei$objf, 1)

    # ... and an ABSOLUTE pin as well, because the relative check alone would
    # pass if a future change degraded the likelihood engine for both methods
    # together.  Measured 19592.53 (foceif) / 19592.63 (focei) at sigdig=3; the
    # window is wide enough for solver-level jitter and far narrower than the
    # +8.63 this test exists to catch.
    expect_equal(.foceif$objf, 19592.5, tolerance = 1e-4)
    expect_equal(.focei$objf, 19592.6, tolerance = 1e-4)

    # the old failure signature was 6 outer evaluations -- one per parameter-ish,
    # independent of the data
    expect_gt(.outerEvals(.foceif), 6)

    # and it should still be doing far less work than the derivative-free default
    expect_lt(.outerEvals(.foceif), .outerEvals(.focei))
  })
})
