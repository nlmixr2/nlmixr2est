# sqrtm(): the matrix square root behind the "|r|" and "|linFim|" repairs.

test_that("sqrtm() is an error for a matrix that is not finite", {
  .msg <- "sqrtm() needs a finite matrix"
  # these gave a 0 x 0 matrix
  expect_error(sqrtm(matrix(c(1, NaN, NaN, 1), 2)), .msg, fixed = TRUE)
  expect_error(sqrtm(matrix(c(1, Inf, Inf, 1), 2)), .msg, fixed = TRUE)
  expect_error(sqrtm(matrix(c(1, 2, NA, 4), 2)), .msg, fixed = TRUE)
  # this one claimed imaginary components, and an infinite diagonal gave an
  # infinite root
  expect_error(sqrtm(matrix(c(NaN, 0, 0, 1), 2)), .msg, fixed = TRUE)
  expect_error(sqrtm(matrix(c(Inf, 0, 0, 1), 2)), .msg, fixed = TRUE)
})

test_that("sqrtm() of a finite matrix is unchanged", {
  .m <- matrix(c(2, 1, 1, 2), 2)
  .r <- sqrtm(.m)
  expect_equal(.r %*% .r, .m, tolerance = 1e-12)
  expect_equal(.r, t(.r))
  # singular: the exact root of a diagonal matrix, the principal root of the
  # rank-one all-ones matrix (J %*% J = 2 J, so J / sqrt(2))
  expect_identical(sqrtm(diag(c(4, 0))), diag(c(2, 0)))
  expect_equal(sqrtm(matrix(1, 2, 2)), matrix(1 / sqrt(2), 2, 2), tolerance = 1e-12)
  expect_identical(dim(sqrtm(matrix(numeric(0), 0, 0))), c(0L, 0L))
  expect_error(sqrtm(matrix(c(1, 2, 2, 1), 2)), "Some components of sqrtm are imaginary.", fixed = TRUE)
})
