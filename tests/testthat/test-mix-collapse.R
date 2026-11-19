test_that(".mixCollapseCheck is silent for well-separated classes", {
  post <- cbind(c(0.99, 0.02, 0.9), c(0.01, 0.98, 0.1))
  ll <- cbind(c(-10, -30, -12), c(-15, -12, -14))
  expect_silent(r <- .mixCollapseCheck(post, ll))
  expect_length(r$empty, 0L)
  expect_equal(nrow(r$same), 0L)
})

test_that(".mixCollapseCheck names an empty class", {
  post <- cbind(c(0.995, 0.999, 0.992), c(0.005, 0.001, 0.008))
  expect_warning(
    r <- .mixCollapseCheck(post),
    "mixture class 2 has no subjects \\(max posterior < 0.01\\); try new starts"
  )
  expect_equal(r$empty, 2L)
  # the boundary itself is not empty
  post[1, ] <- c(0.99, 0.01)
  expect_silent(.mixCollapseCheck(post))
  expect_warning(.mixCollapseCheck(post, emptyTol = 0.0100001), "class 2")
})

test_that(".mixCollapseCheck names an indistinguishable pair", {
  post <- cbind(c(0.9, 0.05, 0.2), c(0.06, 0.03, 0.5), c(0.04, 0.92, 0.3))
  ll <- cbind(c(-10, -30, -12), c(-30, -12, -12.5), c(-30.01, -12.02, -12.49))
  expect_warning(
    r <- .mixCollapseCheck(post, ll),
    "mixture classes 2 and 3 are indistinguishable at the estimates"
  )
  expect_equal(r$same, matrix(c(2L, 3L), ncol = 2L))
  # one subject that tells them apart is enough
  ll[3, 3] <- -12.5 + 0.05
  expect_silent(.mixCollapseCheck(post, ll))
  ll[3, 3] <- -12.5 + 0.0499
  expect_warning(.mixCollapseCheck(post, ll), "classes 2 and 3")
})

test_that(".mixCollapseCheck reports an empty class once, not also as a pair", {
  post <- cbind(c(0.999, 0.998), c(0.001, 0.002))
  ll <- cbind(c(-10, -11), c(-10, -11))
  w <- character(0)
  withCallingHandlers(
    .mixCollapseCheck(post, ll),
    warning = function(e) {
      w <<- c(w, conditionMessage(e))
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(w, "mixture class 2 has no subjects (max posterior < 0.01); try new starts")
})

test_that(".mixCollapseCheck ignores subjects with no evidence for either class", {
  post <- cbind(c(1, 0, 0), c(0, 0.5, 0.5), c(0, 0.5, 0.5))
  # subject 1 belongs to class 1: both other log-posteriors are -Inf
  ll <- log(post) - rep(log(c(1 / 3, 1 / 3, 1 / 3)), each = 3)
  expect_warning(.mixCollapseCheck(post, ll), "classes 2 and 3")
  expect_silent(.mixCollapseCheck(post, NULL))
  expect_silent(.mixCollapseCheck(post[, 1, drop = FALSE]))
})

test_that("mixture warnings stay under 75 characters", {
  expect_lt(nchar("mixture class 10 has no subjects (max posterior < 0.01); try new starts"), 75)
  expect_lt(nchar("mixture classes 10 and 11 are indistinguishable at the estimates"), 75)
})
