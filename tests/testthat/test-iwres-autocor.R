.acDf <- function(id, iwres, time = NULL, ...) {
  if (is.null(time)) {
    time <- stats::ave(seq_along(id), id, FUN = seq_along)
  }
  data.frame(ID = id, TIME = time, IWRES = iwres, ...)
}

test_that("monotone IWRES gives low Durbin-Watson and lag-1 r near 1", {
  .ac <- .iwresAutocor(.acDf(rep(1:3, each = 10), rep(seq(-2, 2, length.out = 10), 3)))
  expect_lt(.ac$durbinWatson, 1.5)
  expect_equal(.ac$lag1R, 1)
  expect_equal(.ac$nSubjects, 3L)
  expect_equal(.ac$nPairs, 27L)
})

test_that("alternating IWRES gives Durbin-Watson above 2 and negative r", {
  .ac <- .iwresAutocor(.acDf(rep(1:4, each = 8), rep(c(1, -1), 16)))
  expect_gt(.ac$durbinWatson, 2)
  expect_equal(.ac$lag1R, -1)
})

test_that("Durbin-Watson matches the within-subject closed form", {
  set.seed(42)
  .id <- rep(1:5, times = c(4, 1, 6, 3, 5))
  .e <- rnorm(length(.id))
  .ac <- .iwresAutocor(.acDf(.id, .e))
  .num <- 0
  .den <- 0
  .x <- .y <- numeric(0)
  for (.i in unique(.id)) {
    .v <- .e[.id == .i]
    if (length(.v) < 2) {
      next
    }
    .num <- .num + sum(diff(.v)^2)
    .den <- .den + sum(.v^2)
    .x <- c(.x, .v[-length(.v)])
    .y <- c(.y, .v[-1])
  }
  expect_equal(.ac$durbinWatson, .num / .den)
  expect_equal(.ac$lag1R, stats::cor(.x, .y))
  expect_equal(.ac$nSubjects, 4L)
  expect_equal(.ac$nPairs, length(.x))
})

test_that("differences never span two subjects", {
  # each subject is constant, so any cross-subject pair would add to the numerator
  .ac <- .iwresAutocor(.acDf(rep(1:3, each = 3), rep(c(-1, 2, 0.5), each = 3)))
  expect_equal(.ac$durbinWatson, 0)
  expect_equal(.ac$nPairs, 6L)
})

test_that("single-observation subjects are skipped", {
  .ac <- .iwresAutocor(.acDf(c(1, 2, 2, 3), c(5, 1, 2, -7)))
  expect_equal(.ac$nSubjects, 1L)
  expect_equal(.ac$nPairs, 1L)
  expect_equal(.ac$durbinWatson, 1 / 5)
  .ac <- .iwresAutocor(.acDf(1:3, c(1, 2, 3)))
  expect_equal(.ac$nPairs, 0L)
  expect_true(is.na(.ac$durbinWatson))
  expect_true(is.na(.ac$lag1R))
})

test_that("NA, dosing and censored rows are skipped", {
  .df <- .acDf(
    rep(1, 6),
    c(1, NA, 2, 100, 3, -100),
    EVID = c(0, 0, 0, 0, 0, 1),
    CENS = c(0, 0, 0, 1, 0, 0)
  )
  .ac <- .iwresAutocor(.df)
  # only 1, 2, 3 remain, paired in time order
  expect_equal(.ac$nPairs, 2L)
  expect_equal(.ac$durbinWatson, 2 / 14)
  # a subject left with a single valid residual contributes nothing
  .df2 <- rbind(.df, .acDf(c(2, 2), c(4, NA), EVID = 0, CENS = 0))
  expect_equal(.iwresAutocor(.df2)$nSubjects, 1L)
})

test_that("rows are ordered by time and endpoints are not paired", {
  .ac <- .iwresAutocor(.acDf(rep(1, 4), c(3, 1, 4, 2), time = c(3, 1, 4, 2)))
  expect_equal(.ac$durbinWatson, 3 / 30)
  # interleaved endpoints: pairs stay within each CMT
  .ac <- .iwresAutocor(.acDf(rep(1, 4), c(1, -5, 1, -5), time = 1:4, CMT = c(1, 2, 1, 2)))
  expect_equal(.ac$durbinWatson, 0)
  expect_equal(.ac$nPairs, 2L)
  expect_equal(.ac$nSubjects, 1L)
})

test_that("a table without IWRES gives NULL", {
  expect_null(.iwresAutocor(data.frame(ID = 1, TIME = 1)))
  expect_null(.iwresAutocor(NULL))
})

test_that("fit$iwresAutocor is computed from the fit table", {
  skip_on_cran()
  .fit <- one.compartment.fit.focei
  .ac <- .fit$iwresAutocor
  expect_named(.ac, c("lag1R", "durbinWatson", "nSubjects", "nPairs"))
  expect_equal(.ac, .iwresAutocor(as.data.frame(.fit)))
  expect_equal(.ac$nSubjects, length(unique(.fit$ID)))
  expect_true(.ac$durbinWatson > 0 && .ac$durbinWatson < 4)
})
