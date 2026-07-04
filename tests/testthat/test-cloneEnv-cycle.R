test_that(".cloneEnv() handles circular environment references without infinite recursion", {
  # Regression test: .cloneEnv() used to recurse forever on circular environment graphs.

  # a -> b -> a  (two-node cycle)
  .a <- new.env(parent = emptyenv())
  .b <- new.env(parent = emptyenv())
  .a$child <- .b
  .b$parent <- .a
  .a$val <- 42L
  .b$val <- 7L

  .cl <- nlmixr2est:::.cloneEnv(.a)

  expect_true(is.environment(.cl))
  expect_false(identical(.cl, .a))          # a real clone, not the original
  expect_equal(.cl$val, 42L)
  expect_true(is.environment(.cl$child))
  expect_false(identical(.cl$child, .b))    # the nested env is cloned too
  expect_equal(.cl$child$val, 7L)
  # cycle is rebuilt pointing at the clone, not duplicated
  expect_identical(.cl$child$parent, .cl)

  # a direct self-reference must also terminate
  .s <- new.env(parent = emptyenv())
  .s$self <- .s
  .s$y <- "z"
  .cs <- nlmixr2est:::.cloneEnv(.s)
  expect_equal(.cs$y, "z")
  expect_identical(.cs$self, .cs)
})
