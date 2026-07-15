## est="advi" thread-count invariance: the per-subject ELBO+gradient core is
## parallelized over subjects (rxControl(cores=)), but a serial id-ordered
## reduction keeps the result bit-for-bit identical to the serial (cores=1) run
## for any thread count.  Uses 2 threads only (CI thread policy).

nmTest({
  one.cmt <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
      d/dt(depot) <- -ka*depot; d/dt(center) <- ka*depot - cl/v*center
      cp <- center/v; cp ~ add(add.sd) })
  }

  .runAdvi <- function(cores, family = "meanField", pointEstimate = TRUE) {
    ctl <- adviControl(iters = 120L, seed = 7L, print = 0L, returnAdvi = TRUE,
                       adviFamily = family, pointEstimate = pointEstimate,
                       rxControl = rxode2::rxControl(cores = cores))
    suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "advi", control = ctl)))
  }

  test_that("mean-field ADVI is bit-for-bit identical for cores=1 vs cores=2", {
    skip_on_cran()
    r1 <- .runAdvi(1L)
    r2 <- .runAdvi(2L)
    expect_identical(r1$theta, r2$theta)
    expect_identical(r1$popOmega, r2$popOmega)
    expect_identical(r1$mu, r2$mu)
    expect_identical(r1$elbo, r2$elbo)
    expect_identical(r1$scale, r2$scale)
  })

  test_that("full-rank ADVI is bit-for-bit identical for cores=1 vs cores=2", {
    skip_on_cran()
    r1 <- .runAdvi(1L, family = "fullRank")
    r2 <- .runAdvi(2L, family = "fullRank")
    expect_identical(r1$theta, r2$theta)
    expect_identical(r1$popOmega, r2$popOmega)
    expect_identical(r1$mu, r2$mu)
    expect_identical(r1$elbo, r2$elbo)
  })

  test_that("full-Bayes ADVI is bit-for-bit identical for cores=1 vs cores=2", {
    skip_on_cran()
    r1 <- .runAdvi(1L, pointEstimate = FALSE)
    r2 <- .runAdvi(2L, pointEstimate = FALSE)
    expect_identical(r1$theta, r2$theta)
    expect_identical(r1$popOmega, r2$popOmega)
    expect_identical(r1$mu, r2$mu)
    expect_identical(r1$elbo, r2$elbo)
    expect_identical(r1$adviCov, r2$adviCov)
  })
})
