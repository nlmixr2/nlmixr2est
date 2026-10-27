test_that(".foceiDropMalformedCov() keeps a proper covariance and one with whole NA rows", {
  .nm <- c("tka", "tcl")
  .pd <- matrix(c(0.04, 0.01, 0.01, 0.09), 2, dimnames = list(.nm, .nm))
  .part <- matrix(c(1, NA, NA, NA), 2, dimnames = list(.nm, .nm))
  for (.ok in list(.pd, .part)) {
    .e <- new.env(parent = emptyenv())
    .e$cov <- .ok
    expect_silent(.foceiDropMalformedCov(.e))
    expect_identical(.e$cov, .ok)
    expect_false(exists("covDebug", envir = .e, inherits = FALSE))
  }
  # no covariance, nothing to check
  expect_silent(.foceiDropMalformedCov(new.env(parent = emptyenv())))
})

test_that(".foceiDropMalformedCov() drops a malformed covariance, says so, and keeps it in covDebug", {
  .nm <- c("tka", "tcl")
  .bad <- matrix(c(1, NA, NA, 1), 2, dimnames = list(.nm, .nm))
  .noNames <- matrix(c(1, 0, 0, 1), 2)
  for (.m in list(.bad, .noNames, "unexpected")) {
    .e <- new.env(parent = emptyenv())
    .e$cov <- .m
    expect_warning(
      suppressMessages(.foceiDropMalformedCov(.e)),
      "^covariance not in proper form, can access value in [$]covDebug$"
    )
    expect_null(.e$cov)
    expect_identical(.e$covDebug, .m)
  }
})

test_that(".foceiDropMalformedCov() has nothing to add to a covariance step that failed", {
  # the nlm family hands over "failed" as the covariance when its step said so
  .e <- new.env(parent = emptyenv())
  .e$cov <- "failed"
  expect_silent(.foceiDropMalformedCov(.e))
  expect_null(.e$cov)
  expect_false(exists("covDebug", envir = .e, inherits = FALSE))
})
