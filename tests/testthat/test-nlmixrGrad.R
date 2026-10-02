test_that("a Gill83 search that ends without an accepted interval restores theta", {
  # x[2] is flat to roundoff, so its search runs out of iterations ("Constant
  # Grad").  It is searched first and df[1] depends on x[2], so an x[2] left at
  # its last probe shows up as a wrong derivative for x[1].
  f <- function(x) (x[1] - 1)^2 + 1 + 1e-9 * x[2]^2 + 1e4 * (x[1] - 1) * x[2]
  g <- nlmixr2Gill83(f, c(1, 2))
  expect_equal(as.character(g$info), c("Good", "Constant Grad"))
  expect_equal(g$df[1], 2e4, tolerance = 1e-7)
})
