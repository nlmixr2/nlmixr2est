# npdeCalc()'s imputation of censored observations, on hand-built normal
# simulations whose conditional truncated distribution is known exactly
nmTest({
  # rxode2's per-thread engines exist once something has been solved
  invisible(rxode2::rxSolve(rxode2::rxode2("y = t"), rxode2::et(1), returnType = "data.frame"))

  # K replicates of nid subjects with n observations each, drawn from N(0, S)
  .npdeSims <- function(S, nid, K) {
    .n <- nrow(S)
    .c <- chol(S)
    .sim <- do.call(
      rbind,
      lapply(seq_len(K), function(k) {
        .z <- matrix(rnorm(nid * .n), nid) %*% .c
        data.frame(simId = k, id = rep(seq_len(nid), each = .n), sim = as.vector(t(.z)))
      })
    )
    .sim$rxLambda <- 1
    .sim$rxYj <- 2
    .sim$rxLow <- 0
    .sim$rxHi <- 1
    .sim
  }

  .npdeImpute <- function(sim, dv, cens, n, limit = NULL) {
    .r <- .Call(
      `_nlmixr2est_npdeCalc`,
      sim,
      dv,
      rep(0L, length(dv)),
      cens,
      limit,
      list(censMethod = 3L, seed = 7L)
    )
    matrix(as.numeric(.r[[1]]), ncol = n, byrow = TRUE)
  }

  test_that("censored values follow the truncated conditional normal", {
    set.seed(1)
    .S <- matrix(c(1, 0.9, 0.9, 0.9, 1, 0.99, 0.9, 0.99, 1), 3)
    .nid <- 1000
    .sim <- .npdeSims(.S, .nid, 300)
    # the first observation measured at 0.5, the others censored below 0
    .imp <- .npdeImpute(.sim, rep(c(0.5, 0, 0), .nid), rep(c(0L, 1L, 1L), .nid), 3)
    expect_equal(.imp[, 1], rep(0.5, .nid))
    expect_true(all(.imp[, 2:3] <= 0))
    # exact, by rejection
    .b <- .S[2:3, 1] / .S[1, 1]
    .v <- .S[2:3, 2:3] - .b %o% .S[1, 2:3]
    .x <- matrix(rnorm(2e6), ncol = 2) %*% chol(.v)
    .x <- sweep(.x, 2, .b * 0.5, "+")
    .x <- .x[.x[, 1] <= 0 & .x[, 2] <= 0, ]
    expect_gt(suppressWarnings(ks.test(.imp[, 2], .x[, 1])$p.value), 0.01)
    expect_gt(suppressWarnings(ks.test(.imp[, 3], .x[, 2])$p.value), 0.01)
    expect_equal(cor(.imp[, 2], .imp[, 3]), cor(.x[, 1], .x[, 2]), tolerance = 0.1)
  })

  test_that("a strongly correlated censored block is not underdispersed", {
    set.seed(2)
    .S <- matrix(0.9999, 3, 3)
    diag(.S) <- 1
    .nid <- 2000
    .sim <- .npdeSims(.S, .nid, 300)
    # censored far above the predictions: effectively untruncated
    .imp <- .npdeImpute(.sim, rep(10, 3 * .nid), rep(1L, 3 * .nid), 3)
    expect_equal(apply(.imp, 2, sd), rep(1, 3), tolerance = 0.06)
    expect_gt(min(cor(.imp)), 0.99)
  })

  test_that("Gibbs leaves a censoring corner it would start on", {
    set.seed(4)
    # eigenvectors with mixed signs, so no single direction is free at the corner
    .u <- cbind(c(1, 1, 1, 1), c(1, -1, 1, -1), c(1, 1, -1, -1), c(1, -1, -1, 1)) / 2
    .S <- .u %*% diag(c(2.5, 1.5, 1, 0.5)) %*% t(.u)
    .nid <- 1000
    .sim <- .npdeSims(.S, .nid, 300)
    # rows 1-3 left censored at 0, row 4 right censored at 0
    .imp <- .npdeImpute(.sim, rep(0, 4 * .nid), rep(c(1L, 1L, 1L, -1L), .nid), 4)
    expect_true(all(.imp[, 1:3] <= 0))
    expect_true(all(.imp[, 4] >= 0))
    .x <- matrix(rnorm(4e6), ncol = 4) %*% chol(.S)
    .x <- .x[.x[, 1] <= 0 & .x[, 2] <= 0 & .x[, 3] <= 0 & .x[, 4] >= 0, ]
    for (.j in 1:4) {
      expect_gt(suppressWarnings(ks.test(.imp[, .j], .x[, .j])$p.value), 0.01, label = .j)
    }
  })

  test_that("a narrow interval does not pin its censored neighbour", {
    set.seed(5)
    .S <- matrix(c(1, 0.5, 0.5, 1), 2)
    .nid <- 2000
    .sim <- .npdeSims(.S, .nid, 300)
    # row 1 interval censored in (0, 1e-6); row 2 effectively untruncated
    .imp <- .npdeImpute(
      .sim,
      rep(c(1e-6, 10), .nid),
      rep(1L, 2 * .nid),
      2,
      limit = rep(c(0, -Inf), .nid)
    )
    expect_true(all(.imp[, 1] >= 0 & .imp[, 1] <= 1e-6))
    # y2 | y1 = 0 is N(0, 0.75)
    expect_equal(sd(.imp[, 2]), sqrt(0.75), tolerance = 0.06)
    expect_equal(mean(.imp[, 2]), 0, tolerance = 0.06)
  })

  test_that("a singular censored block starts inside its support", {
    set.seed(7)
    .nid <- 300
    .sim <- .npdeSims(diag(3), .nid, 300)
    # the third observation is the sum of the first two in every replicate
    .m <- matrix(.sim$sim, ncol = 3, byrow = TRUE)
    .m[, 3] <- .m[, 1] + .m[, 2]
    .sim$sim <- as.vector(t(.m))
    # rows 1-2 left censored at -1, row 3 right censored at -10
    .imp <- suppressWarnings(.npdeImpute(
      .sim,
      rep(c(-1, -1, -10), .nid),
      rep(c(1L, 1L, -1L), .nid),
      3
    ))
    expect_true(all(.imp[, 1:2] <= -1))
    expect_true(all(.imp[, 3] >= -10))
    # strictly below the LOQ, not stuck on it
    expect_lt(max(.imp[, 1:2]), -1 + 1e-6)
    expect_gt(sd(.imp[, 1]), 0.1)
  })

  test_that("a constant simulated row is imputed as its value", {
    set.seed(6)
    .nid <- 500
    .sim <- .npdeSims(diag(2), .nid, 300)
    .sim$sim[c(TRUE, FALSE)] <- 0
    # row 1 (always 0) censored below 1, row 2 censored below 0
    # (the decorrelation itself falls back to eigenvalues and says so)
    .imp <- suppressWarnings(.npdeImpute(.sim, rep(c(1, 0), .nid), rep(1L, 2 * .nid), 2))
    expect_equal(.imp[, 1], rep(0, .nid))
    expect_true(all(.imp[, 2] <= 0))
    # half-normal mean
    expect_equal(mean(.imp[, 2]), -sqrt(2 / pi), tolerance = 0.1)
  })

  test_that("a singular covariance conditions on the measured value", {
    set.seed(3)
    .nid <- 50
    .sim <- .npdeSims(matrix(1, 1, 1), .nid, 300)
    # the second observation duplicates the first in every replicate
    .sim <- .sim[rep(seq_len(nrow(.sim)), each = 2), ]
    # (the decorrelation itself falls back to eigenvalues and says so)
    .imp <- suppressWarnings(.npdeImpute(.sim, rep(c(1, 2), .nid), rep(c(0L, 1L), .nid), 2))
    expect_equal(.imp[, 2], rep(1, .nid), tolerance = 1e-8)
  })
})
