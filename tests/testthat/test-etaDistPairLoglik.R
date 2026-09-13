## The JOINT copula objective for a declared PAIR.
##
## Why it exists: under a copula the marginal parameters' MLE from the joint
## likelihood is NOT the two separate marginal MLEs -- they coincide only at
## rho == 0.  Fitting the marginals separately and carrying a correlation
## alongside them answers a different question than the model asks.
##
## And it makes rho estimable BY THE COPULA DENSITY.  Today the correlation is
## estimated from second moments of the latents, which is a method-of-moments
## estimator, not the MLE; one variant of it (a product-moment correlation of
## the COMBINED latents) is recorded as biased with a fixed point at whatever
## rho it was handed.

.edPairFam <- function() .etaDistFamilyCode("dgamma(2, 2)")

## thetas are (log mean, log relative variance) per marginal, which is how the
## Bauer arms are written: shape = 1/rv, rate = 1/(rv*mean)
.edPairE1 <- c("1/exp(lrv1)", "1/(exp(lrv1)*exp(lm1))")
.edPairE2 <- c("1/exp(lrv2)", "1/(exp(lrv2)*exp(lm2))")
.edPairVars <- c("lm1", "lrv1", "lm2", "lrv2")

.edPairSim <- function(n = 40, rho = 0.5, seed = 1) {
  set.seed(seed)
  .z1 <- stats::rnorm(n)
  .z2 <- rho*.z1 + sqrt(1 - rho^2)*stats::rnorm(n)
  list(e1 = stats::qgamma(stats::pnorm(.z1), 2, 2/5.1),
       e2 = stats::qgamma(stats::pnorm(.z2), 3, 3/50),
       z1 = .z1, z2 = .z2)
}

## an independent R implementation of the same quantity -- deliberately written
## from the definition rather than from the C++, so agreement means something
.edPairRef <- function(e1, e2, rho) {
  .u1 <- pmin(pmax(stats::pgamma(e1, 2, 2/5.1), 1e-15), 1 - 1e-15)
  .u2 <- pmin(pmax(stats::pgamma(e2, 3, 3/50), 1e-15), 1 - 1e-15)
  .q1 <- stats::qnorm(.u1); .q2 <- stats::qnorm(.u2)
  .lc <- -0.5*log(1 - rho^2) -
    (rho^2*(.q1^2 + .q2^2) - 2*rho*.q1*.q2)/(2*(1 - rho^2))
  sum(stats::dgamma(e1, 2, 2/5.1, log = TRUE) +
        stats::dgamma(e2, 3, 3/50, log = TRUE) + .lc)
}

.edPairObj <- function(theta, e1, e2, rho, rhoIdx = -1L, vars = .edPairVars) {
  rxEtaDistPairLoglikTest_(.edPairFam(), .edPairFam(),
                           .edPairE1, .edPairE2, vars, theta,
                           matrix(0, length(e1), 0), e1, e2,
                           rho, rhoIdx, rep(1, length(e1)))
}

test_that("the joint objective matches an independent implementation", {
  .s <- .edPairSim()
  .th <- c(log(5.1), log(1/2), log(50), log(1/3))
  for (.rho in c(0, 0.3, 0.5, -0.4)) {
    .o <- .edPairObj(.th, .s$e1, .s$e2, .rho)
    expect_equal(length(.o), 1L)
    expect_equal(.o[1], .edPairRef(.s$e1, .s$e2, .rho), tolerance = 1e-10)
  }
  ## rho == 0 must reduce to the two marginals with no copula term at all
  expect_equal(.edPairObj(.th, .s$e1, .s$e2, 0)[1],
               sum(stats::dgamma(.s$e1, 2, 2/5.1, log = TRUE)) +
                 sum(stats::dgamma(.s$e2, 3, 3/50, log = TRUE)),
               tolerance = 1e-10)
})

test_that("the gradient agrees with central differences", {
  ## The pair gradient differences the ASSEMBLED density rather than using the
  ## analytic d(log p)/da times a differenced da/dtheta, because the copula term
  ## reaches the arguments through z = qnorm(F(x; a)) and no analytic dF/da is
  ## available.  So this check is the only thing standing behind it.
  .s <- .edPairSim()
  .th <- c(log(5.1), log(1/2), log(50), log(1/3))
  .a <- attr(.edPairObj(.th, .s$e1, .s$e2, 0.5), "grad")
  expect_false(is.null(.a))
  for (.i in seq_along(.th)) {
    .h <- 1e-5*max(abs(.th[.i]), 1)
    .tp <- .th; .tp[.i] <- .tp[.i] + .h
    .tm <- .th; .tm[.i] <- .tm[.i] - .h
    .n <- (.edPairObj(.tp, .s$e1, .s$e2, 0.5)[1] -
             .edPairObj(.tm, .s$e1, .s$e2, 0.5)[1])/(2*.h)
    expect_equal(.a[.i], .n, tolerance = 1e-5, info = .edPairVars[.i])
  }
})

test_that("rho is estimable by the copula density, and unbiasedly", {
  ## rhoIdx >= 0 moves the correlation into the parameter vector on the atanh
  ## scale.  First: the two forms must agree when they describe the same rho.
  .s <- .edPairSim()
  .th <- c(log(5.1), log(1/2), log(50), log(1/3))
  .v2 <- c(.edPairVars, "rxRho")
  expect_equal(.edPairObj(c(.th, atanh(0.5)), .s$e1, .s$e2, 0, 4L, .v2)[1],
               .edPairObj(.th, .s$e1, .s$e2, 0.5)[1], tolerance = 1e-12)
  ## its gradient too
  .g <- attr(.edPairObj(c(.th, atanh(0.5)), .s$e1, .s$e2, 0, 4L, .v2), "grad")
  .h <- 1e-5
  .f <- function(.a) .edPairObj(c(.th, .a), .s$e1, .s$e2, 0, 4L, .v2)[1]
  expect_equal(.g[5], (.f(atanh(0.5) + .h) - .f(atanh(0.5) - .h))/(2*.h),
               tolerance = 1e-5)
  ## and the maximizer tracks the SAMPLE, which is what an MLE does -- this
  ## particular draw has cor(z1,z2) = 0.627, not the 0.5 it was simulated at.
  ##
  ## NEAR the sample correlation, not equal to it: with the marginal parameters
  ## at their true values the reconstructed z = qnorm(F(eta)) IS the simulated
  ## latent, so the only gap left is that the Gaussian-copula MLE of rho solves
  ## a cubic rather than being the sample correlation.  Measured 0.660 against
  ## 0.627 at n = 40.  Absolute, because a relative tolerance on a correlation
  ## says something different at 0.05 than at 0.9.
  .mle <- function(e1, e2) {
    tanh(stats::optimize(function(.a)
      -.edPairObj(c(.th, .a), e1, e2, 0, 4L, .v2)[1], c(-3, 3))$minimum)
  }
  expect_lt(abs(.mle(.s$e1, .s$e2) - stats::cor(.s$z1, .s$z2)), 0.06)
  ## and it is on the same side of 0.5 as the sample is, i.e. it followed the
  ## data rather than sitting where it started
  expect_gt(.mle(.s$e1, .s$e2), 0.55)
  ## unbiased over replicates, and shrinking like sqrt(n).  Measured at 200
  ## replicates: n=40 gave 0.5131 (sd 0.0956), n=400 gave 0.4964 (sd 0.0345).
  skip_on_cran()
  .est <- vapply(seq_len(60), function(.i) {
    .r <- .edPairSim(n = 200, rho = 0.5, seed = 1000 + .i)
    .mle(.r$e1, .r$e2)
  }, numeric(1))
  expect_equal(mean(.est), 0.5, tolerance = 0.03)
})

test_that("the objective declines rather than guessing", {
  .s <- .edPairSim()
  .th <- c(log(5.1), log(1/2), log(50), log(1/3))
  ## an expression outside the grammar -> zero-length, not a wrong number
  expect_equal(length(rxEtaDistPairLoglikTest_(
    .edPairFam(), .edPairFam(), c("sin(lrv1)", "1"), .edPairE2,
    .edPairVars, .th, matrix(0, 40, 0), .s$e1, .s$e2, 0.5, -1L, rep(1, 40))), 0L)
  ## an eta outside the family's support -> declines, does not return a partial
  ## sum over the records that happened to work
  .bad <- .s$e1; .bad[3] <- -1
  expect_equal(length(.edPairObj(.th, .bad, .s$e2, 0.5)), 0L)
  ## a degenerate correlation has no density
  expect_equal(length(.edPairObj(.th, .s$e1, .s$e2, 1)), 0L)
})
