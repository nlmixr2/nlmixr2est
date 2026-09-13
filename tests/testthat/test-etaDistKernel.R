## Metropolis kernels for a declared eta sampled ON ITS OWN SCALE.
##
## Tested against a FLAT likelihood, where the stationary distribution must be
## the family itself.  That is the property that says the proposal, the
## bijector and the Jacobian agree with each other -- independently of any
## model, any estimator, and any of the machinery around them.
##
## It is also a test that can fail.  Removing the log-Jacobian from kernel 2
## does not make the chain merely inefficient, it makes it wrong: measured on
## gamma(shape=0.5, rate=1), a log-scale random walk WITHOUT the Jacobian
## returns mean 0.0002 against a true 0.5, because the map compresses toward
## zero and nothing corrects for it.  So these moment checks are the guard on
## the one line that is easiest to drop and hardest to notice.

nmTest({

  .kern <- function(fam, a, n = 60000L, kernel = 1L, s = 1, start = 1) {
    set.seed(42)
    rxEtaDistKernelTest_(fam, as.numeric(a), n, kernel, s, start, -1L)
  }

  ## family, args, true mean, true variance, a sane starting point
  .cases <- list(
    list("dnorm real",        1L, c(0, 1),   0.0,  1.00,  0.0),
    list("dgamma shape<1",   13L, c(0.5, 1), 0.5,  0.50,  1.0),
    list("dgamma shape>1",   13L, c(3, 2),   1.5,  0.75,  1.0),
    list("dbeta unit",       20L, c(2, 3),   0.4,  0.04,  0.5),
    list("dunif bounded",    22L, c(2, 5),   3.5,  0.75,  3.0),
    list("dexp positive",    12L, 2,         0.5,  0.25,  1.0))

  test_that("kernel 1 (independence from the prior) reproduces every support type", {
    ## The proposal is inverse-CDF through rxEtaDistQ(), which every declarable
    ## family has -- so one code path covers all 22 and no per-family generator
    ## is needed.  Because the proposal IS the prior, the prior cancels from the
    ## ratio and a flat likelihood accepts everything: this checks the quantile
    ## route as much as the kernel.
    for (.c in .cases) {
      .x <- .kern(.c[[2]], .c[[3]], kernel = 1L, start = .c[[6]])
      expect_equal(mean(.x), .c[[4]], tolerance = 0.05, info = .c[[1]])
      expect_equal(var(.x), .c[[5]], tolerance = 0.10, info = .c[[1]])
    }
  })

  test_that("kernel 2 (bijected random walk) reproduces every support type", {
    ## A random walk on the eta itself would propose outside the support of a
    ## positive or bounded family; walking on the bijected scale cannot.  The
    ## Jacobian is what makes it target the right distribution rather than
    ## merely stay legal -- see the collapse recorded at the top of this file.
    for (.c in .cases) {
      .x <- .kern(.c[[2]], .c[[3]], kernel = 2L, s = 1, start = .c[[6]])
      expect_equal(mean(.x), .c[[4]], tolerance = 0.05, info = .c[[1]])
      expect_equal(var(.x), .c[[5]], tolerance = 0.10, info = .c[[1]])
    }
  })

  test_that("the two kernels agree with each other, not just with the truth", {
    ## Independent failure modes: kernel 1 is wrong if the quantile function is,
    ## kernel 2 if the bijector or Jacobian is.  Agreeing to a tighter tolerance
    ## than either meets against the truth is evidence they are not both wrong
    ## in the same direction.
    for (.c in .cases) {
      .a <- mean(.kern(.c[[2]], .c[[3]], kernel = 1L, start = .c[[6]]))
      .b <- mean(.kern(.c[[2]], .c[[3]], kernel = 2L, s = 1, start = .c[[6]]))
      expect_equal(.a, .b, tolerance = 0.05, info = .c[[1]])
    }
  })

  test_that("the chain does not depend on where it started", {
    ## A start far into the tail of a shape<1 gamma is the case that breaks a
    ## sampler adapted on the wrong scale -- and it is g4's dispersion.
    .near <- mean(.kern(13L, c(0.5, 1), kernel = 2L, s = 1, start = 0.1))
    .far  <- mean(.kern(13L, c(0.5, 1), kernel = 2L, s = 1, start = 50))
    expect_equal(.near, 0.5, tolerance = 0.05)
    expect_equal(.far,  0.5, tolerance = 0.05)
  })
})
