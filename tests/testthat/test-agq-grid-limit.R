# The nAGQ^neta grid is checked before expand.grid() builds it (#1412)
test_that(".agq() refuses an oversized grid without building it", {
  expect_error(.agq(neta = 12, nAGQ = 5), "244,140,625 nodes")
  expect_error(.agq(neta = 40, nAGQ = 9), "exceeds the limit")
})

test_that(".agq() small grid is unchanged and silent", {
  expect_silent(g <- .agq(neta = 2, nAGQ = 3))
  expect_equal(g$n, 9L)
  expect_equal(dim(g$x), c(9L, 2L))
  expect_equal(dim(g$w), c(9L, 2L))
  expect_true(g$first)
  expect_equal(sum(apply(g$w, 1, prod)), 1)
})

test_that(".agq() warns once above the soft threshold", {
  w <- character(0)
  g <- withCallingHandlers(
    .agq(neta = 6, nAGQ = 5),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(w, 1L)
  expect_match(w, "15,625 nodes")
  expect_lt(nchar(w), 75L)
  expect_equal(g$n, 15625L)
  expect_silent(.agq(neta = 6, nAGQ = 5, warnNodes = Inf))
})
