## The declared eta's CDF, and the JOINT density of a declared pair on the eta
## scale under a Gaussian copula.
##
## These are what let the direct route carry a CORRELATED declared block.  A
## Gaussian copula over non-normal marginals is exactly `eta = Q(phi(z))`, so
## refusing such a block was refusing the only thing that gives the correlation
## meaning; written on the eta scale it is just another prior term.
##
## It also moves the inverse CDF out of the model.  Decoding runs once per
## OBSERVATION; evaluating this runs once per eta per SUBJECT.

nmTest({

  .fams <- c("dnorm","stdNormal","studentT","dcauchy","doubleExponential","dlogis",
             "gumbel","dlnorm","dchisq","invChiSquare","scaledInvChiSquare","dexp",
             "dgamma","invGamma","dweibull","frechet","rayleigh","pareto",
             "paretoType2","dbeta","betaProportion","dunif")
  .args <- list(c(0,1), numeric(0), c(5,0,1), c(0,1), c(0,1), c(0,1), c(0,1),
                c(0,1), 3, 3, c(3,1), 1, c(0.5,1), c(2,1), c(2,1), c(2,1), 1,
                c(1,2), c(0,1,2), c(2,3), c(0.5,5), c(0,1))

  test_that("the CDF is the exact inverse of the quantile function", {
    ## P(Q(u)) == u, family by family.  This is the test that catches a
    ## REVERSED TAIL -- invChiSquare, scaledInvChiSquare, invGamma and pareto
    ## all define Q through the upper tail, so a naive CDF is the complement of
    ## what the quantile says.  That mismatch is silent: it does not error, it
    ## produces a wrong copula.
    for (.k in seq_len(22L)) {
      for (.u in c(0.01, 0.1, 0.5, 0.9, 0.99)) {
        .r <- rxEtaDistPTest_(.k, .u, as.numeric(.args[[.k]]))
        expect_equal(unname(.r[["pBack"]]), .u, tolerance = 1e-9,
                     info = paste(.fams[.k], "at u =", .u))
      }
    }
  })

  test_that("the copula's marginal z survives the far tail of a shape<1 gamma", {
    ## qnorm(F(eta)) is the composition this route rests on, and gamma with
    ## shape 1/rv = 0.5 -- g4's dispersion, the arm with no interior mode -- is
    ## where it is most likely to saturate.  It must come back as qnorm(u).
    for (.u in c(0.001, 0.1, 0.5, 0.9, 0.999)) {
      .r <- rxEtaDistPTest_(13L, .u, c(0.5, 1))
      expect_equal(unname(.r[["z"]]), stats::qnorm(.u), tolerance = 1e-6,
                   info = paste("u =", .u))
    }
  })

  test_that("the joint density integrates to one", {
    ## The check that a density is a density.  A mis-signed exponent or a
    ## dropped -0.5*log(1-rho^2) still returns finite, plausible numbers and
    ## still has a maximum in the right place -- it just does not integrate to
    ## 1, which is the only thing that catches it.
    .a <- c(0.5, 1)                       # gamma shape 0.5
    .n <- 200L
    .u <- (seq_len(.n) - 0.5)/.n
    .q <- stats::qgamma(.u, .a[1], .a[2])
    .lm <- stats::dgamma(.q, .a[1], .a[2], log = TRUE)
    for (.rho in c(0, 0.5, -0.7)) {
      .tot <- 0
      for (.i in seq_len(.n)) {
        for (.j in seq_len(.n)) {
          .lp <- rxEtaDistPairLogDTest_(13L, .q[.i], .a, 13L, .q[.j], .a, .rho)
          .tot <- .tot + exp(.lp - .lm[.i] - .lm[.j])
        }
      }
      expect_equal(.tot/.n^2, 1, tolerance = 2e-3,
                   info = paste("rho =", .rho))
    }
  })

  test_that("rho = 0 is EXACTLY the sum of the marginals", {
    ## Independence must cost nothing.  If the copula term leaks in at rho = 0
    ## every uncorrelated declared model silently gets a wrong prior.
    .a <- c(0.5, 1)
    expect_equal(rxEtaDistPairLogDTest_(13L, 0.8, .a, 13L, 1.2, .a, 0),
                 stats::dgamma(0.8, .a[1], .a[2], log = TRUE) +
                   stats::dgamma(1.2, .a[1], .a[2], log = TRUE),
                 tolerance = 1e-12)
  })

  test_that("a degenerate correlation has no density rather than a wrong one", {
    .a <- c(2, 1)
    expect_true(is.infinite(rxEtaDistPairLogDTest_(13L, 1, .a, 13L, 1, .a, 1)))
    expect_true(is.infinite(rxEtaDistPairLogDTest_(13L, 1, .a, 13L, 1, .a, -1)))
  })

  test_that("the joint works across DIFFERENT families, not just a matched pair", {
    ## The marginals are independent choices -- a gamma clearance and a
    ## lognormal volume is an ordinary model -- so the pair must not assume
    ## they match.
    .g <- c(0.5, 1); .l <- c(0, 1)
    .v <- rxEtaDistPairLogDTest_(13L, 0.8, .g, 8L, 1.3, .l, 0.4)
    expect_true(is.finite(.v))
    ## and at rho = 0 it is still the sum of the two different marginals
    expect_equal(rxEtaDistPairLogDTest_(13L, 0.8, .g, 8L, 1.3, .l, 0),
                 stats::dgamma(0.8, .g[1], .g[2], log = TRUE) +
                   stats::dlnorm(1.3, .l[1], .l[2], log = TRUE),
                 tolerance = 1e-12)
  })
})
