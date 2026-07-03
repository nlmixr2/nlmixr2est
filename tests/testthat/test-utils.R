test_that(".collectWarn", {
  expect_equal(
    .collectWarn(
      {
        warning("A")
        1
      }, lst = TRUE),
    list(1, warning = "A", error = NULL)
  )
  # warnings are not duplicated
  expect_equal(
    .collectWarn(
      {
        warning("A")
        warning("A")
        1
      }, lst = TRUE),
    list(1, warning = "A", error = NULL)
  )
  # Non-lst raises the warning
  expect_warning(
    check <- .collectWarn(
      {
        warning("A")
        warning("A")
        1
      }),
    regexp = "A"
  )
  expect_equal(check, 1)
  # Errors are not captured by default
  expect_error(
    .collectWarn(
      {
        warning("A")
        stop("B")
        1
      }),
    regexp = "B"
  )
  # Errors are not captured by default
  expect_error(
    check <- .collectWarn(
      {
        warning("A")
        stop("B")
        1
      }, lst = TRUE, collectErr = TRUE),
    regexp = NA
  )
  expect_equal(check, list(NULL, warning = "A", error = "B"))
  # Inner try() blocks still work — when expr ultimately succeeds, the
  # captured "internal" errors are discarded so they cannot leak into
  # the caller's error list (regression: previously every error fired
  # inside expr was returned, breaking SAEM's covariance recovery).
  expect_equal(
    .collectWarn(
      {
        try(stop("internal"), silent = TRUE)
        42
      }, lst = TRUE, collectErr = TRUE),
    list(42, warning = NULL, error = NULL)
  )
  # When the expression fails after an on.exit handler raises a second
  # error (issue 607 pattern), every message observed along the error
  # chain is preserved.  on.exit must be inside a function frame for
  # this to trigger; that mirrors the real-world rxUiGet.foceiHdEta
  # case in R/focei.R.
  .foceiHdEtaLike <- function() {
    on.exit(stop("Aborted calculation"))
    stop("none of the predictions depend on 'ETA'")
  }
  expect_error(
    check <- .collectWarn(.foceiHdEtaLike(), lst = TRUE, collectErr = TRUE),
    regexp = NA
  )
  expect_true(is.null(check[[1]]))
  expect_true("Aborted calculation" %in% check$error)
  expect_true("none of the predictions depend on 'ETA'" %in% check$error)
})
