## Inner-block weight-injection hook: nlmixrSetInnerWeightFn registers a callback
## that FOCEI invokes right after a subject's etas are written to par_ptr and
## before its inner solve (on every eta-set, including finite-difference
## perturbations).  A plugin (nlmixr2nn) uses it to overwrite the subject's
## par_ptr weight block with individual weights W_i = f(lW, etaW); FOCEI's FD
## eta-sensitivity then captures d(f)/d(etaW) through the injected weights.
## Here we validate the mechanism fires with the current etas during a fit.

test_that("inner-weight injection hook fires per subject eta-set during a fit", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  .old <- rxode2::getRxThreads(); on.exit(rxode2::setRxThreads(.old), add = TRUE)
  rxode2::setRxThreads(1L)

  d <- nlmixr2data::theo_sd
  m <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7; eta.cl ~ 0.1 })
    model({
      ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  .Call("_nlmixr2est_registerTestInnerWt", PACKAGE = "nlmixr2est")
  on.exit(.Call("_nlmixr2est_removeTestInnerWt", PACKAGE = "nlmixr2est"), add = TRUE)

  f <- suppressWarnings(suppressMessages(
    nlmixr2(m, d, est = "focei",
            control = foceiControl(print = 0L, maxOuterIterations = 2L,
                                   maxInnerIterations = 3L, calcTables = FALSE))))

  res <- .Call("_nlmixr2est_getTestInnerWt", PACKAGE = "nlmixr2est")
  expect_gt(res[1], 0)                 # hook fired during the inner problem
  expect_true(is.finite(res[2]))       # etas passed through
})
