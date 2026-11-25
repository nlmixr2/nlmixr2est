test_that("a declared argument with unary minus is evaluated as R evaluates it", {
  # `-a - b` was read as -(a - b): dgamma(rate = exp(-lclrv - lclm)) became a
  # gamma with mean exp(-lclm), and the direct route estimated lclm as minus
  # its value
  .f <- rxEtaDistLoglikTest_
  .eta <- c(2, 3, 5, 7, 4.5)
  .th <- c(0.7, 1.6)
  .ll <- function(exprs) {
    .v <- .f(13L, exprs, c("lclrv", "lclm"), .th, matrix(0, 5L, 0L), .eta, rep(1, 5))
    if (length(.v) == 0L) NA_real_ else .v[1]
  }
  .r <- sum(stats::dgamma(.eta, shape = exp(-0.7), rate = exp(-0.7 - 1.6), log = TRUE))
  expect_equal(.ll(c("exp(-lclrv)", "exp(-lclrv - lclm)")), .r, tolerance = 1e-10)
  expect_equal(.ll(c("exp(-lclrv)", "exp(-lclrv-lclm)")), .r, tolerance = 1e-10)
  expect_equal(.ll(c("exp(-lclrv)", "exp(-(lclrv+lclm))")), .r, tolerance = 1e-10)
  expect_equal(.ll(c("exp(-lclrv)", "exp(-lclrv)/exp(lclm)")), .r, tolerance = 1e-10)
  # -x^2 is -(x^2), and unary minus binds tighter than * and /
  .r2 <- sum(stats::dgamma(.eta, shape = exp(-(0.7^2)), rate = exp(-0.7 * 1.6), log = TRUE))
  expect_equal(.ll(c("exp(-lclrv^2)", "exp(-lclrv*lclm)")), .r2, tolerance = 1e-10)
})
