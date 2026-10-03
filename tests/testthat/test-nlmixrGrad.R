test_that("a Gill83 search that ends without an accepted interval restores theta", {
  # x[2] is flat to roundoff, so its search runs out of iterations ("Constant
  # Grad").  It is searched first and df[1] depends on x[2], so an x[2] left at
  # its last probe shows up as a wrong derivative for x[1].
  f <- function(x) (x[1] - 1)^2 + 1 + 1e-9 * x[2]^2 + 1e4 * (x[1] - 1) * x[2]
  g <- nlmixr2Gill83(f, c(1, 2))
  expect_equal(as.character(g$info), c("Good", "Constant Grad"))
  expect_equal(g$df[1], 2e4, tolerance = 1e-7)
  # with the last parameter left out, the base objective was never evaluated
  g <- nlmixr2Gill83(f, c(1, 2), which = c(TRUE, FALSE))
  expect_equal(as.character(g$info), c("Good", "Not Assessed"))
  expect_equal(g$f[1], 1)
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

test_that("nlmixr2Gill83() and nlmixr2Hess() use the gill* arguments", {
  g <- nlmixr2Gill83(
    sin,
    1,
    gillRtol = 1e-6,
    gillK = 4L,
    gillStep = 3,
    gillFtol = 1e-3
  )
  expect_equal(g$gillRtol, 1e-6)
  expect_equal(g$gillK, 4L)
  expect_equal(g$gillStep, 3)
  expect_equal(g$gillFtol, 1e-3)
  # an accepted interval is 2*sqrt(|f|*gillRtol/|f''|)
  expect_equal(as.character(g$info), "Good")
  expect_equal(g$hf, 2 * sqrt(abs(sin(1)) * 1e-6 / abs(g$df2)))
  expect_false(identical(nlmixr2Hess(1, sin, gillRtol = 1e-4), nlmixr2Hess(1, sin)))
  expect_error(nlmixr2Gill83(sin, 1, gillK = -1L), "gillK")
})

test_that("gillK = 0 determines no step size", {
  # the search would start from gillStep * 2 * (1 + |x|) * sqrt(epsA / (1 + |f|)),
  # epsA = |f| * gillRtol; that is what is reported, unsearched
  .f <- sin(1)
  .epsA <- abs(.f) * sqrt(.Machine$double.eps)
  .h0 <- 2 * (2 * (1 + 1) * sqrt(.epsA / (1 + abs(.f))))
  g <- nlmixr2Gill83(sin, 1, gillK = 0L)
  expect_identical(as.character(g$info), "Not Assessed")
  expect_identical(g$hf, .h0)
  expect_identical(g$hphi, .h0)
  expect_true(is.na(g$df) && is.na(g$df2) && is.na(g$err))
  # a function flat in x is not searched either
  g <- nlmixr2Gill83(function(x) 1 + 1e-9 * x^2, 2, gillK = 0L)
  expect_identical(as.character(g$info), "Not Assessed")
  # nlmixr2Hess differences with that interval: |x| * rEpsC + aEpsC = .h0
  .h <- 1 * (.h0 * 0.5) + .h0 * 0.5
  .H <- (-sin(1 + 2 * .h) + 16 * sin(1 + .h) - 30 * .f + 16 * sin(1 - .h) - sin(1 - 2 * .h)) /
    (12 * .h * .h)
  expect_equal(nlmixr2Hess(1, sin, gillK = 0L), matrix(.H))
  # and so does the first nlmixr2GradFun gradient, a forward difference
  gf <- nlmixr2GradFun(sin, gillK = 0L, print = 0)
  gf$eval(1)
  expect_equal(gf$grad(1), (sin(1 + .h) - .f) / .h)
})

test_that("nlmixr2GradFun() gradients leave the point alone", {
  gf <- nlmixr2GradFun(
    function(x) if (x[1] > 1) NA_real_ else sum(x^2) + 1,
    print = 0
  )
  x <- c(1, 1)
  gf$eval(x)
  gf$grad(x) # the first gradient is the Gill search
  # x[1]'s forward leg is NA, so it takes the backward difference
  g1 <- gf$grad(x)
  g2 <- gf$grad(x)
  expect_identical(x, c(1, 1))
  expect_identical(g1, g2)
  expect_equal(g1[2], 2, tolerance = 1e-3)
})
