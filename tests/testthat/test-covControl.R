# The finite-difference covariance options foceiControl() and rsControl() share.

test_that("foceiControl() and rsControl() take the same finite-difference covariance options", {
  # gillStepCov is the factor the Gill search grows its step by
  expect_error(foceiControl(gillStepCov = 0.5), "gillStepCov")
  expect_error(rsControl(gillStepCov = 0.5), "gillStepCov")
  expect_identical(foceiControl(gillStepCov = 1)$gillStepCov, 1)
  expect_identical(rsControl(gillStepCov = 1)$gillStepCov, 1)
  # covSmall is a single finite threshold
  expect_error(foceiControl(covSmall = Inf), "covSmall")
  expect_error(rsControl(covSmall = Inf), "covSmall")
  expect_error(foceiControl(covSmall = c(1e-5, 1e-6)), "covSmall")
  expect_error(rsControl(covSmall = c(1e-5, 1e-6)), "covSmall")
  # the flags are logical or 0/1, which foceiControl() stores
  expect_identical(
    unclass(rsControl(covGillF = 0L, rmatNorm = 1, smatNorm = TRUE)),
    list(covGillF = 0L, rmatNorm = 1, smatNorm = TRUE)
  )
  expect_identical(
    foceiControl(covGillF = 0L, rmatNorm = 1, smatNorm = TRUE)[c("covGillF", "rmatNorm", "smatNorm")],
    list(covGillF = 0L, rmatNorm = 1L, smatNorm = 1L)
  )
  expect_error(rsControl(rmatNorm = 2), "rmatNorm")
  expect_error(foceiControl(rmatNorm = 2), "rmatNorm")
})
