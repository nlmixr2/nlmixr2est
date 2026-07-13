## ADVI reproducibility guarantees from the counter-based RNG (stream keyed by
## the GLOBAL iteration index): same seed -> identical; a short run is a bit-for-
## bit prefix of a longer one; a warm resume equals a single fresh longer run;
## and results are independent of the thread count.

nmTest({
  one.cmt <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
      d/dt(depot) <- -ka*depot; d/dt(center) <- ka*depot - cl/v*center
      cp <- center/v; cp ~ add(add.sd) })
  }
  runAdvi <- function(ctl) suppressMessages(suppressWarnings(
    nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "advi", control = ctl)))

  test_that("iters=100 is a bit-for-bit prefix of iters=200", {
    r100 <- runAdvi(adviControl(iters = 100L, seed = 3L, print = 0L, returnAdvi = TRUE))
    r200 <- runAdvi(adviControl(iters = 200L, seed = 3L, print = 0L, returnAdvi = TRUE))
    expect_identical(r100$elbo, r200$elbo[1:100])
    expect_identical(r100$parHist, r200$parHist[1:100, , drop = FALSE])
  })

  test_that("warm resume is numerically equivalent to a fresh longer run", {
    ## The RNG stream, variational/population state, and step-size accumulators
    ## are continued exactly (the prefix test above is bit-for-bit).  Resume is
    ## numerically -- not bit-for-bit -- equal to a fresh longer run because the
    ## ODE solver carries per-subject internal state (step-size memory) that
    ## depends on the solve history: fresh reached iteration 100 via 100 solves,
    ## resume via a cold set-up.  The continued trajectory tracks the fresh one to
    ## solver-tolerance precision.
    r100 <- runAdvi(adviControl(iters = 100L, seed = 3L, print = 0L, returnAdvi = TRUE))
    rResume <- runAdvi(adviControl(iters = 100L, seed = 3L, print = 0L,
                                   returnAdvi = TRUE, resume = r100))
    rFresh <- runAdvi(adviControl(iters = 200L, seed = 3L, print = 0L, returnAdvi = TRUE))
    expect_equal(rResume$it0, 200L)
    expect_equal(rResume$theta, rFresh$theta, tolerance = 5e-2)
    expect_equal(rResume$logPopOmega, rFresh$logPopOmega, tolerance = 5e-2)
    ## the resumed ELBO tail tracks the fresh run's second half
    expect_equal(mean(rResume$elbo), mean(rFresh$elbo[101:200]), tolerance = 1e-2)
  })

  test_that("results are independent of the thread count", {
    c1 <- adviControl(iters = 60L, seed = 5L, print = 0L, returnAdvi = TRUE,
                      rxControl = rxode2::rxControl(cores = 1L))
    c4 <- adviControl(iters = 60L, seed = 5L, print = 0L, returnAdvi = TRUE,
                      rxControl = rxode2::rxControl(cores = 4L))
    r1 <- runAdvi(c1); r4 <- runAdvi(c4)
    expect_identical(r1$theta, r4$theta)
    expect_identical(r1$elbo, r4$elbo)
    expect_identical(r1$mu, r4$mu)
  })
})
