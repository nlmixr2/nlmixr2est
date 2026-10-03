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

test_that("every Gill step factor is a finite number of at least 1", {
  # gillStep (the fit's gradient search), gillStepCov and gillStepCovLlik (the
  # covariance step) all grow the Gill (1983) step by multiplying by the factor
  # and shrink it by dividing
  for (.n in c("gillStep", "gillStepCov", "gillStepCovLlik")) {
    .below <- paste0("Assertion on '", .n, "' failed: Element 1 is not >= 1.")
    .inf <- paste0("Assertion on '", .n, "' failed: Must be finite.")
    expect_error(do.call(foceiControl, setNames(list(0.5), .n)), .below, fixed = TRUE)
    expect_error(do.call(foceiControl, setNames(list(0), .n)), .below, fixed = TRUE)
    expect_error(do.call(foceiControl, setNames(list(Inf), .n)), .inf, fixed = TRUE)
    expect_identical(do.call(foceiControl, setNames(list(1), .n))[[.n]], 1)
  }
  expect_error(rsControl(gillStepCov = Inf), "Assertion on 'gillStepCov' failed: Must be finite.", fixed = TRUE)
  # nlmixr2Gill83() makes the checks foceiControl() makes
  .f <- function(x) sum(x^2)
  expect_error(
    nlmixr2Gill83(.f, c(1, 2), gillStep = 0.5),
    "Assertion on 'gillStep' failed: Element 1 is not >= 1.",
    fixed = TRUE
  )
  expect_error(
    nlmixr2Gill83(.f, c(1, 2), gillStep = Inf),
    "Assertion on 'gillStep' failed: Must be finite.",
    fixed = TRUE
  )
})
