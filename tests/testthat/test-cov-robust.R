test_that(".nlmixr2RobustCov handles full-rank matrices like the plain QR solve", {
  set.seed(42)
  X <- matrix(rnorm(40), ncol = 4)
  expected <- solve(crossprod(X))
  got <- nlmixr2est:::.nlmixr2RobustCov(X)
  expect_equal(unname(got), unname(expected), tolerance = 1e-6)
})

test_that(".nlmixr2RobustCov degrades only the ill-identified parameter(s)", {
  set.seed(42)
  X <- matrix(rnorm(40), ncol = 4)
  # Zero out one column -- unidentifiable, should not blow up the others
  X[, 3] <- 0
  got <- nlmixr2est:::.nlmixr2RobustCov(X)
  expect_equal(dim(got), c(4L, 4L))
  expect_true(all(is.na(got[3, ])))
  expect_true(all(is.na(got[, 3])))
  # The well-identified 3x3 submatrix should match a direct solve
  keep <- c(1, 2, 4)
  expected <- solve(crossprod(X[, keep, drop = FALSE]))
  expect_equal(unname(got[keep, keep]), unname(expected), tolerance = 1e-6)
})

test_that(".nlmixr2RobustCov returns all-NA when nothing is identifiable", {
  X <- matrix(0, nrow = 10, ncol = 3)
  got <- nlmixr2est:::.nlmixr2RobustCov(X)
  expect_true(all(is.na(got)))
  expect_equal(dim(got), c(3L, 3L))
})

test_that(".nlmixr2CholPartial validates only the finite submatrix", {
  covm <- diag(c(1, 2, 3))
  dimnames(covm) <- list(c("a", "b", "c"), c("a", "b", "c"))
  # Fully finite, positive-definite: behaves like chol()
  expect_false(inherits(nlmixr2est:::.nlmixr2CholPartial(covm), "try-error"))

  covm2 <- covm
  covm2[2, ] <- NA_real_
  covm2[, 2] <- NA_real_
  # NA row/col for one parameter, rest still positive-definite: still valid
  expect_false(inherits(nlmixr2est:::.nlmixr2CholPartial(covm2), "try-error"))

  covm3 <- covm
  covm3[1, 1] <- -1 # not positive-definite even ignoring NA
  expect_true(inherits(nlmixr2est:::.nlmixr2CholPartial(covm3), "try-error"))

  covm4 <- matrix(NA_real_, 2, 2)
  expect_true(inherits(nlmixr2est:::.nlmixr2CholPartial(covm4), "try-error"))
})
