## The declared eta on ITS OWN scale: bounds, bijector, log-Jacobian, and the
## two derivatives of the log density with respect to the ETA.
##
## Everything else in the etaDist C++ works on the family's ARGUMENTS.  This is
## the other axis, and it is what a sampler proposing the eta directly needs and
## what an inner MAP needs in order to use a non-Gaussian prior -- the Gaussian
## quadratic `eta' Omega^-1 eta` and its curvature `Omega^-1` are the special
## case of `-2 log p` and `-d2 log p/d(eta)2`.

nmTest({

  .es <- function(fam, x, a) rxEtaDistEtaScaleTest_(fam, x, as.numeric(a))

  test_that("support comes from the ARGUMENTS, not lotri's support column", {
    ## The column says what a family is morally, and for a bijector it is wrong
    ## in three places.  Taking it at face value puts an identity map on a
    ## bounded eta, and every proposal outside the range comes back -Inf with
    ## nothing saying why -- so these three are the reason the bounds are
    ## computed per family rather than read off the table.
    .u <- .es(22L, 3, c(2, 5))          # dunif  -- column says "real"
    expect_equal(unname(.u[["lo"]]), 2)
    expect_equal(unname(.u[["hi"]]), 5)
    .p <- .es(18L, 4, c(3, 2))          # pareto -- column says "positive"
    expect_equal(unname(.p[["lo"]]), 3) # y_min, not 0
    expect_true(is.infinite(.p[["hi"]]))
    .p2 <- .es(19L, 2, c(1, 2, 3))      # paretoType2 -- column says "nonneg"
    expect_equal(unname(.p2[["lo"]]), 1)
    ## and the ordinary positive family is still [0, Inf)
    .g <- .es(13L, 1, c(0.5, 1))
    expect_equal(unname(.g[["lo"]]), 0)
    expect_true(is.infinite(.g[["hi"]]))
  })

  test_that("the bijector round-trips on every declarable family", {
    ## x -> u -> x must be the identity for all 22, because a proposal is made
    ## in u and scored in x; any family where that drifts silently biases the
    ## chain rather than erroring.
    .a <- list(c(0,1), numeric(0), c(5,0,1), c(0,1), c(0,1), c(0,1), c(0,1),
               c(0,1), 3, 3, c(3,1), 1, c(0.5,1), c(2,1), c(2,1), c(2,1), 1,
               c(1,2), c(0,1,2), c(2,3), c(0.5,5), c(0,1))
    .x <- c(0.3,0.3,0.3,0.3,0.3,0.3,0.3, 1.2,2.0,0.5,0.5, 0.7,0.8,1.5,1.1,
            1.3,0.9,2.5,1.0, 0.4,0.4,0.5)
    for (.k in seq_len(22L)) {
      .r <- .es(.k, .x[.k], .a[[.k]])
      expect_equal(unname(.r[["xBack"]]), .x[.k], tolerance = 1e-8,
                   info = paste("family code", .k))
      expect_true(is.finite(.r[["u"]]), info = paste("family code", .k))
      expect_true(is.finite(.r[["logJac"]]), info = paste("family code", .k))
    }
  })

  test_that("the log-Jacobian matches the map it belongs to", {
    ## unbounded: identity, so 0
    expect_equal(unname(.es(1L, 0.3, c(0,1))[["logJac"]]), 0)
    ## lower-bounded: x = lo + exp(u), so log|dx/du| = u
    .g <- .es(13L, 0.8, c(0.5, 1))
    expect_equal(unname(.g[["logJac"]]), unname(.g[["u"]]))
    ## bounded both sides: finite-difference the inverse map through the round trip
    .lo <- 0; .hi <- 1; .x <- 0.4
    .b <- .es(20L, .x, c(2, 3))
    .u <- unname(.b[["u"]]); .h <- 1e-6
    .xu <- function(u) .lo + (.hi - .lo)/(1 + exp(-u))
    expect_equal(unname(.b[["logJac"]]),
                 log((.xu(.u + .h) - .xu(.u - .h))/(2*.h)), tolerance = 1e-6)
  })

  test_that("the eta derivatives match closed forms", {
    ## dnorm(0,1): d = -x, d2 = -1
    .r <- .es(1L, 0.3, c(0, 1))
    expect_equal(unname(.r[["dlogp"]]), -0.3, tolerance = 1e-5)
    expect_equal(unname(.r[["d2logp"]]), -1, tolerance = 1e-5)
    ## dexp(rate): d = -rate, d2 = 0
    .r <- .es(12L, 0.7, 2)
    expect_equal(unname(.r[["dlogp"]]), -2, tolerance = 1e-5)
    expect_equal(unname(.r[["d2logp"]]), 0, tolerance = 1e-5)
    ## dgamma(shape, rate): d = (shape-1)/x - rate, d2 = -(shape-1)/x^2
    .r <- .es(13L, 0.8, c(0.5, 1))
    expect_equal(unname(.r[["dlogp"]]), (0.5 - 1)/0.8 - 1, tolerance = 1e-5)
    expect_equal(unname(.r[["d2logp"]]), -(0.5 - 1)/0.8^2, tolerance = 1e-5)
  })

  test_that("a shape < 1 gamma reports POSITIVE curvature, and is not hidden", {
    ## This is the Pumas caution in one number.  At g4's dispersion the shape is
    ## 1/rv = 0.5, the log density is CONVEX, and d2 log p/d(eta)2 is positive --
    ## so the prior contributes NEGATIVE curvature to an inner Hessian and there
    ## is no interior mode to expand a Laplace approximation around.  The
    ## primitive must report that rather than smooth it away, so whatever
    ## consumes it can refuse or fall back instead of silently producing a
    ## non-positive-definite Hessian.
    expect_true(.es(13L, 0.8, c(0.5, 1))[["d2logp"]] > 0)
    ## shape > 1 has an interior mode and the usual negative curvature
    expect_true(.es(13L, 2.0, c(3, 1))[["d2logp"]] < 0)
    ## and the boundary case, exponential, is exactly flat in curvature
    expect_equal(unname(.es(13L, 1.0, c(1, 1))[["d2logp"]]), 0, tolerance = 1e-5)
  })
})
