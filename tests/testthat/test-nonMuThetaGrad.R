nmTest({
  # The exact-gradient core (src/nonMuThetaGrad.cpp) is shared by every
  # estimator that takes this step, so it is checked on its own rather than only
  # through a fit.  A fit needs a model whose sensitivity peer fits the solve
  # pool, which is incidental to whether the derivatives are right -- and while
  # that was unavailable the arithmetic went unverified.

  # the Gaussian branch differentiates exactly phi0NormalSSR's per-observation
  # term (src/saem.cpp):  0.5*((y-f)/g)^2 + log(g),  g = ares + bres*|f|
  .obj <- function(y, f0, a, b, dfdth) {
    function(th) {
      .f <- f0 + sum(dfdth * th)
      .g <- a + b * abs(.f)
      0.5 * ((y - .f) / .g)^2 + log(.g)
    }
  }
  .fd <- function(fn, th, h = 1e-6) {
    vapply(seq_along(th), function(i) {
      .tp <- th; .tm <- th; .tp[i] <- .tp[i] + h; .tm[i] <- .tm[i] - h
      (fn(.tp) - fn(.tm)) / (2 * h)
    }, numeric(1))
  }

  test_that("the Gaussian score matches numerical differentiation", {
    withr::with_seed(1, {
      for (.k in 1:6) {
        .y <- rnorm(1, 5, 2)
        .f0 <- .y + rnorm(1, 0, 0.5)
        .a <- runif(1, 0.05, 0.5)
        # every other case is proportional error, where the residual scale moves
        # with the prediction -- dropping that second term is the classic way to
        # get this quietly wrong on a prop() model, so it must be exercised
        .b <- if (.k %% 2 == 0) runif(1, 0.05, 0.3) else 0
        .d <- rnorm(3)
        .g <- .a + .b * abs(.f0)
        .r <- nlmixr2est:::nonMuGradAccumTest_(0L, .y, .f0, .g, .b * sign(.f0), .d, 1)
        .num <- .fd(.obj(.y, .f0, .a, .b, .d), rep(0, 3))
        expect_equal(as.numeric(.r$score), .num, tolerance = 1e-5)
      }
    })
  })

  test_that("the log-likelihood score is -d(loglik)/d(theta)", {
    # f IS the per-observation loglik and the objective is -sum(f)
    .d <- c(1.5, -2.5, 0.75)
    .r <- nlmixr2est:::nonMuGradAccumTest_(2L, 0, 1.23, 0, 0, .d, 1)
    expect_equal(as.numeric(.r$score), -.d)
    # BHHH information is the outer product of that per-observation score
    expect_equal(as.numeric(.r$info), as.numeric(outer(.d, .d)))
  })

  test_that("the information matrix is symmetric and non-negative definite", {
    .d <- c(0.4, -1.1, 2.0)
    .r <- nlmixr2est:::nonMuGradAccumTest_(0L, 4.0, 3.5, 0.3, 0.1, .d, 1)
    .i <- matrix(as.numeric(.r$info), 3, 3)
    expect_equal(.i, t(.i))
    expect_true(all(eigen(.i, only.values = TRUE)$values >= -1e-10))
  })

  test_that("a non-finite observation contributes nothing", {
    .d <- c(1, 2, 3)
    .r <- nlmixr2est:::nonMuGradAccumTest_(0L, NA_real_, 3.5, 0.3, 0.0, .d, 1)
    expect_equal(as.numeric(.r$score), rep(0, 3))
    # a zero or negative residual SD is not a valid Gaussian either
    .r0 <- nlmixr2est:::nonMuGradAccumTest_(0L, 4.0, 3.5, 0.0, 0.0, .d, 1)
    expect_equal(as.numeric(.r0$score), rep(0, 3))
  })
})
