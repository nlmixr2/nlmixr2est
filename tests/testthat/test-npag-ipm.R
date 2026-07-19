# Burke interior-point solver (npIpmBurke) vs golden fixtures captured from the
# Pmetrics Rust reference (see design/npag/README.md).  Core, always-run: pure
# numerics, no model fit.

.npReadBurkeGolden <- function() {
  .f <- testthat::test_path("fixtures", "npag", "ipm", "burke-golden.csv")
  .lines <- readLines(.f)
  .cases <- list()
  for (.ln in .lines) {
    .p <- strsplit(.ln, ",", fixed = TRUE)[[1]]
    .name <- .p[1]; .kind <- .p[2]
    .vals <- as.numeric(.p[-(1:3)])
    if (is.null(.cases[[.name]])) .cases[[.name]] <- list(psi = list(), weights = NULL, objective = NULL)
    if (.kind == "psi") {
      .cases[[.name]]$psi[[length(.cases[[.name]]$psi) + 1L]] <- .vals
    } else if (.kind == "weights") {
      .cases[[.name]]$weights <- .vals
    } else if (.kind == "objective") {
      .cases[[.name]]$objective <- .vals
    }
  }
  lapply(.cases, function(.c) {
    .c$psi <- do.call(rbind, .c$psi)  # subjects x support points
    .c
  })
}

test_that("npIpmBurke matches the Pmetrics golden fixtures", {
  .cases <- .npReadBurkeGolden()
  expect_true(length(.cases) >= 4L)
  for (.name in names(.cases)) {
    .c <- .cases[[.name]]
    .res <- npIpmBurke(.c$psi)
    # weights: non-negative, sum to 1, and match the golden values
    expect_true(all(.res$weights >= -1e-12), info = paste(.name, "non-negative"))
    expect_equal(sum(.res$weights), 1, tolerance = 1e-8, info = paste(.name, "sums to 1"))
    expect_equal(.res$weights, .c$weights, tolerance = 1e-6, info = paste(.name, "weights"))
    expect_equal(.res$objective, .c$objective, tolerance = 1e-6, info = paste(.name, "objective"))
  }
})

test_that("npIpmBurke recovers uniform weights for symmetric matrices", {
  # identity -> 1/n; all-ones -> 1/ncol (analytic, backend-independent)
  .id <- npIpmBurke(diag(5))
  expect_equal(.id$weights, rep(1 / 5, 5), tolerance = 1e-8)
  .ones <- npIpmBurke(matrix(1, 4, 6))
  expect_equal(.ones$weights, rep(1 / 6, 6), tolerance = 1e-8)
})

test_that("npIpmBurke handles negatives (abs) and errors on non-finite", {
  # negative entries are coerced to absolute value: |−5| dominates
  .m <- matrix(1, 2, 3); .m[1, 1] <- -5
  .w <- npIpmBurke(.m)$weights
  expect_true(.w[1] > .w[2] && .w[1] > .w[3])
  .bad <- matrix(1, 3, 3); .bad[1, 1] <- NA_real_
  expect_error(npIpmBurke(.bad), "finite")
})
