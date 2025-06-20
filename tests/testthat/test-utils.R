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
  expect_equal(check, list(NULL, warning = "A", error = "doWithOneRestart(return(expr), restart): B"))
})
