## Outer-problem NN-training hook: FOCEI hands a registered R callback a
## method-agnostic per-observation snapshot -- [rx_pred_, <ODE states>, <full
## calc_lhs row>] -- on each real objective evaluation.  The state columns carry
## the NN-weight forward-sensitivity states rx_sw when the model has them; the
## lhs columns carry the NN output g, output transforms, rx_drdg, error pieces
## and any covariates emitted as rx_<cov>_.  The method-specific d(LL)/d(f) comes
## from the contribution hook, not this matrix.  Here we validate the mechanism
## (fires, correct shape incl. the lhs block, f column sane) on a plain model.

test_that("FOCEI outer NN hook passes the method-agnostic prediction+state matrix", {
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
  ## now [rx_pred_, states, lhs_vars]: f + >=2 states + the appended calc_lhs row
  expect_gt(cap$dim[2], 3)
  ## rx_pred_ (col 1) is also present within the appended lhs block
  .matchesF <- which(vapply(seq_len(ncol(cap$mat)),
    function(j) isTRUE(all.equal(cap$mat[, j], cap$mat[, 1], tolerance = 1e-8)),
    logical(1)))
  expect_gt(length(.matchesF), 1)
  ## column 1 is the predicted value f: finite, non-negative, and in a plausible
  ## range for the observed data (method-agnostic -- no dv/r sent)
  expect_true(all(is.finite(cap$mat[, 1])) && all(cap$mat[, 1] >= 0))
  expect_lt(max(cap$mat[, 1]), 2 * max(d$DV[d$EVID == 0]))
})
