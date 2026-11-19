# .getCachedFit() (helper-zzz-fits.R) writes the fixture cache that parallel
# workers and concurrent runs of one checkout share: it has to write it whole.

.fitCacheTest <- function(dir, version = "v1") {
  .f <- .getCachedFit
  environment(.f) <- list2env(
    list(.fitCacheDir = dir, .fitOnCran = FALSE, .fitCacheVersion = version),
    parent = environment(.getCachedFit)
  )
  .f
}

test_that(".getCachedFit() writes the cache whole and leaves no partial file", {
  .dir <- withr::local_tempdir()
  .f <- .fitCacheTest(.dir)
  .runs <- new.env(parent = emptyenv())
  .runs$n <- 0L
  .fit <- function() {
    .runs$n <- .runs$n + 1L
    stats::lm(dist ~ speed, datasets::cars)
  }
  .first <- suppressMessages(.f("lm", .fit, "lm.rds"))
  expect_identical(list.files(.dir, all.files = TRUE, no.. = TRUE), "lm.rds")
  .cached <- readRDS(file.path(.dir, "lm.rds"))
  expect_identical(.cached$version, "v1")
  expect_equal(stats::coef(.cached$fit), stats::coef(.first))
  # a second call reads the file instead of fitting again
  .second <- suppressMessages(.f("lm", .fit, "lm.rds"))
  expect_identical(.runs$n, 1L)
  expect_equal(stats::coef(.second), stats::coef(.first))
  # another version refits and replaces it
  .third <- suppressMessages(.fitCacheTest(.dir, "v2")("lm", .fit, "lm.rds"))
  expect_identical(.runs$n, 2L)
  expect_identical(readRDS(file.path(.dir, "lm.rds"))$version, "v2")
  expect_identical(list.files(.dir, all.files = TRUE, no.. = TRUE), "lm.rds")
})

test_that(".getCachedFit() leaves nothing behind when the cache cannot be written", {
  .dir <- withr::local_tempdir()
  # a directory where the file belongs: it cannot be read or replaced
  dir.create(file.path(.dir, "lm.rds"))
  .w <- character()
  .fit <- NULL
  .w <- testthat::capture_warnings(
    .fit <- suppressMessages(
      .fitCacheTest(.dir)("lm", function() stats::lm(dist ~ speed, datasets::cars), "lm.rds")
    )
  )
  expect_s3_class(.fit, "lm")
  expect_true(any(grepl("Failed to save cache for lm", .w, fixed = TRUE)))
  expect_identical(list.files(.dir, all.files = TRUE, no.. = TRUE), "lm.rds")
  expect_true(dir.exists(file.path(.dir, "lm.rds")))
})
