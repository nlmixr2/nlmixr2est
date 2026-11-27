test_that("burn-in trend test: false-trend rate under AR(1) noise", {
  .ar1 <- function(n, rho) {
    .x <- numeric(n)
    .e <- rnorm(n)
    .x[1] <- .e[1] / sqrt(1 - rho^2)
    for (.i in 2:n) {
      .x[.i] <- rho * .x[.i - 1] + .e[.i]
    }
    .x
  }
  .rate <- function(rho, nRep = 2000L) {
    .m <- rxode2::rxWithSeed(1401, vapply(seq_len(nRep), function(i) .ar1(100L, rho), numeric(100)))
    mean(saemBurnTrendTest(.m, 5L)$trend)
  }
  .r0 <- .rate(0)
  expect_gt(.r0, 0.03)
  expect_lt(.r0, 0.07)
  # batch means of 20 keep heavy autocorrelation's inflation modest
  .r9 <- .rate(0.9)
  expect_gt(.r9, .r0)
  expect_lt(.r9, 0.15)
})

test_that("burn-in trend test: ramp, constant and non-finite columns", {
  .n <- 100L
  .noise <- rxode2::rxWithSeed(1, rnorm(.n, sd = 0.1))
  .x <- cbind(
    ramp = seq(0, 1, length.out = .n) + .noise,
    flat = .noise,
    constant = rep(3, .n),
    negInf = rep(-Inf, .n),
    hasNa = c(rep(1, .n - 1L), NA),
    hasInf = c(.noise[-1], Inf),
    line = seq(0, 1, length.out = .n)
  )
  .r <- saemBurnTrendTest(.x, 5L)
  expect_equal(.r$trend, c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
  expect_equal(.r$t[3:4], c(0, 0))
  expect_true(all(is.finite(.r$t[1:2])))
  # the slope is tested on the trailing nBatch * floor(n / nBatch) rows
  .y <- c(rep(100, 3), .x[, "flat"])
  expect_false(saemBurnTrendTest(matrix(.y), 5L)$trend)
})
