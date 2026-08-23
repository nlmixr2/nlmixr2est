# Phase 3b.3: the substituted carry gradient, checked numerically against
# finite differences on eta through a real rxSolve() of the generated inner
# model.  The generated text itself is checked in
# test-focei-lincmt-carry-emit.R; shared fixtures in helper-lincmt-carry.R.
# Every test here needs an rxode2 with the carry sentinels and skips cleanly
# on a released rxode2 without them.

.carrySolveInner <- function(m, pars, ev) {
  function(...) {
    .p <- pars
    .set <- list(...)
    for (.n in names(.set)) .p[.n] <- .set[[.n]]
    rxode2::rxSolve(m, params = .p, events = ev, returnType = "data.frame")
  }
}

.carryRelErr <- function(got, fd) {
  max(abs(got - fd) / (abs(fd) + 1e-8))
}

test_that("substituted gradient matches FD-on-eta; naive build shows the bug", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  .ev <- .carryEv()
  .pars <- c(`THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = 0.5,
             `ETA[1]` = 0.3)
  .sC <- .carrySetControl(.carryUiCov(), "auto")$foceiEnv
  .sN <- .carrySetControl(.carryUiCov(), "none")$foceiEnv
  .mC <- suppressWarnings(rxode2::rxode2(.sC$..inner))
  .mN <- suppressWarnings(rxode2::rxode2(.sN$..inner))
  .solveC <- .carrySolveInner(.mC, .pars, .ev)
  .solveN <- .carrySolveInner(.mN, .pars, .ev)
  .h <- 1e-5
  .fd <- (.solveC(`ETA[1]` = 0.3 + .h)$rx_pred_ -
            .solveC(`ETA[1]` = 0.3 - .h)$rx_pred_) / (2 * .h)
  .r0 <- .solveC()
  .rn <- .solveN()
  # predictions themselves agree between the two builds
  expect_equal(.r0$rx_pred_, .rn$rx_pred_, tolerance = 1e-10)
  expect_lt(.carryRelErr(.r0$rx__sens_rx_pred__BY_ETA_1___, .fd), 1e-6)
  # the naive gradient is measurably wrong once wt changes within the subject
  expect_gt(.carryRelErr(.rn$rx__sens_rx_pred__BY_ETA_1___, .fd), 1e-3)
})

test_that("two simultaneous pairs (cl and v slots) both match FD", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  .ev <- .carryEv()
  .pars <- c(`THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = 0.5,
             `ETA[1]` = 0.3, `ETA[2]` = -0.2)
  .m <- suppressWarnings(rxode2::rxode2(.carryUiTwoPair()$foceiEnv$..inner))
  .solve <- .carrySolveInner(.m, .pars, .ev)
  .h <- 1e-5
  .r0 <- .solve()
  .fd1 <- (.solve(`ETA[1]` = 0.3 + .h)$rx_pred_ -
             .solve(`ETA[1]` = 0.3 - .h)$rx_pred_) / (2 * .h)
  .fd2 <- (.solve(`ETA[2]` = -0.2 + .h)$rx_pred_ -
             .solve(`ETA[2]` = -0.2 - .h)$rx_pred_) / (2 * .h)
  expect_lt(.carryRelErr(.r0$rx__sens_rx_pred__BY_ETA_1___, .fd1), 1e-6)
  # eta.v needs the row-local direct term on top of the amounts carry
  expect_lt(.carryRelErr(.r0$rx__sens_rx_pred__BY_ETA_2___, .fd2), 1e-6)
})

test_that("proportional error chains d(rx_r_)/d(eta) through the carried sensitivity", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  .s <- nlmixr2est::nlmixr2(.carryModProp)$foceiEnv
  .inner <- .s$..inner
  expect_true(grepl("rx__sens_rx_r__BY_ETA_1___=.*rx__sens_rx_pred__BY_ETA_1___",
                    .inner))
  .ev <- .carryEv()
  .pars <- c(`THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = 0.1,
             `ETA[1]` = 0.3)
  .m <- suppressWarnings(rxode2::rxode2(.inner))
  .solve <- .carrySolveInner(.m, .pars, .ev)
  .h <- 1e-5
  .r0 <- .solve()
  .fdR <- (.solve(`ETA[1]` = 0.3 + .h)$rx_r_ -
             .solve(`ETA[1]` = 0.3 - .h)$rx_r_) / (2 * .h)
  expect_lt(.carryRelErr(.r0$rx__sens_rx_r__BY_ETA_1___, .fdR), 1e-6)
})
