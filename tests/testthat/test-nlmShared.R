# The nlm-family covariance from the optimizer's Hessian (R/nlmShared.R).

test_that("a positive-definite Hessian is inverted as is", {
  .h <- matrix(c(4, 1, 1, 3), 2)
  .r <- .nlmCovFromHessian(.h)
  expect_identical(.r$type, "r")
  expect_identical(.r$r, .h)
  expect_null(.r$warning)
})

test_that("an indefinite Hessian is repaired as |r|, else as the nearest positive-definite matrix", {
  .h <- matrix(c(1, 2, 2, 1), 2) # eigenvalues 3, -1
  .r <- .nlmCovFromHessian(.h)
  expect_identical(.r$type, "|r|")
  expect_equal(.r$r, sqrtm(.h %*% .h))
  expect_equal(eigen(.r$r, symmetric = TRUE, only.values = TRUE)$values, c(3, 1))
  expect_identical(.r$warning, "R matrix is not positive definite; corrected as \"|r|\"")
})

test_that("a numerically singular Hessian is not repaired", {
  # |R| and the nearest positive-definite matrix would both invert rounding
  # noise: nmNearPD() floors the zero eigenvalue to 2e-8 (a variance of 5e7)
  .sing <- "R matrix is singular; covariance step failed"
  for (.h in list(
    diag(c(2, 0, -1)),
    diag(c(2, -1e-17, 1)), # |R| passes as positive definite
    tcrossprod(1:4) # rank one; |R| fails and R+ was installed
  )) {
    .r <- .nlmCovFromHessian(.h)
    expect_identical(.r$type, "failed")
    expect_null(.r$r)
    expect_identical(.r$warning, .sing)
  }
  # a well-conditioned indefinite Hessian is still repaired
  expect_identical(.nlmCovFromHessian(diag(c(2, -1, 1)))$type, "|r|")
})

test_that("a Hessian that cannot be repaired gives no covariance", {
  .r <- .nlmCovFromHessian(matrix(c(1, NA, NA, 1), 2))
  expect_identical(.r$type, "failed")
  expect_null(.r$r)
  expect_identical(.r$warning, "R matrix is not finite; covariance step failed")
  .r <- .nlmCovFromHessian(matrix(c(1, Inf, Inf, 1), 2))
  expect_identical(.r$warning, "R matrix is not finite; covariance step failed")
  # the zero-filled Hessian of a failed solve
  .r <- .nlmCovFromHessian(matrix(0, 2, 2))
  expect_identical(.r$type, "failed")
  expect_identical(.r$warning, "R matrix is singular; covariance step failed")
})

nmTest({
  .pk <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(centr) <- ka * depot - cl / v * centr
      cp <- centr / v
      cp ~ add(add.sd)
    })
  }

  test_that("a derivative-free fit with an indefinite Hessian gets a repaired positive-definite covariance", {
    .fit <- suppressMessages(nlmixr2(.pk, nlmixr2data::theo_sd, est = "bobyqa", control = bobyqaControl(print = 0L)))
    .h <- .fit$env$bobyqa$r
    expect_lt(min(eigen(.h, symmetric = TRUE, only.values = TRUE)$values), 0)
    expect_identical(.fit$covMethod, "|r|")
    expect_gt(min(eigen(.fit$cov, symmetric = TRUE, only.values = TRUE)$values), 0)
    expect_equal(unname(.fit$env$bobyqa$cov.scaled), unname(solve(sqrtm(.h %*% .h))), tolerance = 1e-8)
    expect_true("R matrix is not positive definite; corrected as \"|r|\"" %in% .fit$runInfo)
    expect_null(.fit$env$bobyqa$covWarning)
  })
})
