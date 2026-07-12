## External likelihood-contribution API: a contributor registers begin/obs/end
## hooks that nlmixr2est cycles in series inside likInner0's per-observation
## loop.  Validated with a test contributor (src/inner.cpp) that (a) counts the
## begin/obs/end calls and (b) adds a constant c to each observation's log-
## likelihood -- which must shift the objective by exactly -2*c*nObs, since a
## constant LL addition does not move the inner optimum.

test_that("likelihood-contribution hooks fire per observation and fold into the objective", {
  skip_on_cran()

  one.compartment <- function() {
    ini({
      tka <- fix(0.45); tcl <- fix(1); tv <- fix(3.45)
      eta.ka ~ fix(0.6); eta.cl ~ fix(0.3); eta.v ~ fix(0.1)
      add.sd <- fix(0.7)
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  .old <- rxode2::getRxThreads(); on.exit(rxode2::setRxThreads(.old), add = TRUE)
  rxode2::setRxThreads(1L)   # test contributor uses global accumulators

  .nObs <- sum(theo_sd$EVID == 0)
  .nsub <- length(unique(theo_sd$ID))

  f0 <- .nlmixr(one.compartment, theo_sd, est = "focei", control = foceiControl(print = 0L))

  cc <- 0.01
  .Call("_nlmixr2est_registerTestContrib", PACKAGE = "nlmixr2est")
  on.exit(.Call("_nlmixr2est_removeTestContrib", PACKAGE = "nlmixr2est"), add = TRUE)
  .Call("_nlmixr2est_setTestContribAddLL", cc, PACKAGE = "nlmixr2est")
  f1 <- .nlmixr(one.compartment, theo_sd, est = "focei", control = foceiControl(print = 0L))
  res <- .Call("_nlmixr2est_getTestContrib", PACKAGE = "nlmixr2est")
  names(res) <- c("nObs", "sumDLLdf", "sumErr", "sumF", "nBegin", "nEnd",
                  "nlhs", "lhsHasF")

  ## a constant c added per obs shifts the objective by exactly -2*c*nObs
  expect_equal(f1$objf - f0$objf, -2 * cc * .nObs, tolerance = 1e-3)

  ## the obs hook fires once per observation within each subject bracket
  expect_gt(res[["nObs"]], 0)
  expect_equal(res[["nObs"]] / res[["nBegin"]], .nObs / .nsub)  # obs per subject
  ## begin/end are balanced
  expect_equal(res[["nBegin"]], res[["nEnd"]])
  ## the solved LHS row is exposed to the obs hook and holds this obs's values
  ## (the prediction f appears among the lhs entries at every observation)
  expect_gt(res[["nlhs"]], 0)
  expect_equal(res[["lhsHasF"]], res[["nObs"]])

  ## after removing the contributor, a fit does not invoke it
  .Call("_nlmixr2est_removeTestContrib", PACKAGE = "nlmixr2est")
  .Call("_nlmixr2est_registerTestContrib", PACKAGE = "nlmixr2est")  # resets counters
  .Call("_nlmixr2est_removeTestContrib", PACKAGE = "nlmixr2est")
  f2 <- .nlmixr(one.compartment, theo_sd, est = "focei", control = foceiControl(print = 0L))
  res2 <- .Call("_nlmixr2est_getTestContrib", PACKAGE = "nlmixr2est")
  expect_equal(res2[[1]], 0)   # nObs: not invoked
})
