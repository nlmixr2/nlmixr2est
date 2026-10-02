test_that("a Gill83 search that ends without an accepted interval restores theta", {
  # x[2] is flat to roundoff, so its search runs out of iterations ("Constant
  # Grad").  It is searched first and df[1] depends on x[2], so an x[2] left at
  # its last probe shows up as a wrong derivative for x[1].
  f <- function(x) (x[1] - 1)^2 + 1 + 1e-9 * x[2]^2 + 1e4 * (x[1] - 1) * x[2]
  g <- nlmixr2Gill83(f, c(1, 2))
  expect_equal(as.character(g$info), c("Good", "Constant Grad"))
  expect_equal(g$df[1], 2e4, tolerance = 1e-7)
})

test_that("nlmixr2Gill83() and nlmixr2Hess() never modify the caller's vector", {
  p <- c(1, 2)
  p0 <- p + 0
  h <- nlmixr2Hess(p, function(x) (x[1] - 1)^2 + 1)
  expect_identical(p, p0)
  expect_equal(h, matrix(c(2, 0, 0, 0), 2), tolerance = 1e-6)
  # an error part-way through a search used to leave the vector at the probe
  f <- function(x) if (x[1] > 1.0001) stop("outside") else sum(x^2)
  expect_error(nlmixr2Gill83(f, p), "outside")
  expect_identical(p, p0)
  # every call gets its own vector: one the objective keeps is not changed by
  # the probes after it, so only the base evaluation was taken at p
  acc <- new.env(parent = emptyenv())
  acc$x <- list()
  f2 <- function(x) {
    acc$x[[length(acc$x) + 1L]] <- x
    sum((x - c(1, 2))^2) + 1
  }
  p2 <- c(1, 2)
  h <- nlmixr2Hess(p2, f2)
  expect_identical(p2, p0)
  expect_equal(h, diag(2, 2), tolerance = 1e-6)
  expect_identical(sum(vapply(acc$x, identical, logical(1), p0)), 1L)
})
