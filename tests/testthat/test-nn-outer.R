## Outer-problem NN-training hook: FOCEI hands a registered R callback a matrix
## with one row per observation -- [id, time, dv, pred, r, <every solved state>]
## -- on each real objective evaluation.  The state columns carry the NN-weight
## forward-sensitivity states (rx_sw) when the model has them; here we validate
## the mechanism (fires, correct shape, dv/pred sane) on a plain model.

test_that("FOCEI outer NN hook passes the per-observation prediction+state matrix", {
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

  cap <- new.env(); cap$n <- 0L
  fn <- function(mat) { cap$n <- cap$n + 1L; cap$dim <- dim(mat); cap$mat <- mat }
  .Call("_nlmixr2est_setNnOuterFn", fn, PACKAGE = "nlmixr2est")
  on.exit(.Call("_nlmixr2est_setNnOuterFn", NULL, PACKAGE = "nlmixr2est"), add = TRUE)

  f <- suppressWarnings(suppressMessages(
    nlmixr2(m, d, est = "focei",
            control = foceiControl(print = 0L, maxOuterIterations = 2L,
                                   maxInnerIterations = 3L, calcTables = FALSE))))

  .nObs <- sum(d$EVID == 0)
  expect_gt(cap$n, 0)                                   # hook fired on real objective evals
  expect_equal(cap$dim[1], .nObs)                       # one row per observation
  expect_gt(cap$dim[2], 5)                              # 5 fixed cols + >=1 state
  ## dv column (4th, 1-based) is the observed data
  expect_equal(sort(cap$mat[, 3]), sort(d$DV[d$EVID == 0]), tolerance = 1e-6)
  ## prediction column (5th) is finite and non-negative
  expect_true(all(is.finite(cap$mat[, 4])) && all(cap$mat[, 4] >= 0))
})
