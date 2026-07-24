## est="advi" end-to-end (mean-field, point-estimate): the raw C++ loop (ELBO
## trend, reproducibility) and the finalized nlmixr2FitData (objective, tables,
## covariance) assembled via nlmixr2CreateOutputFromUi at the ADVI estimates.

nmTest({
  one.cmt <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
      d/dt(depot) <- -ka*depot; d/dt(center) <- ka*depot - cl/v*center
      cp <- center/v; cp ~ add(add.sd) })
  }

  test_that("est='advi' raw loop: ELBO increases and is reproducible", {
    ctl <- adviControl(iters = 150L, seed = 7L, print = 0L, returnAdvi = TRUE)
    res <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "advi", control = ctl)))
    expect_s3_class(res, "nlmixr2advi")
    e <- res$elbo; n <- length(e); d <- max(1L, n %/% 10L)
    expect_gt(mean(e[(n - d + 1L):n]), mean(e[1:d]))   # ELBO trend up
    expect_true(all(is.finite(res$theta)))
    expect_true(all(res$popOmega > 0))

    ## same seed -> identical
    res2 <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "advi", control = ctl)))
    expect_identical(res$theta, res2$theta)
    expect_identical(res$elbo, res2$elbo)
    expect_identical(res$mu, res2$mu)
  })

  test_that("est='advi' assembles a full nlmixr2FitData (objf + tables + cov)", {
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "advi",
              control = adviControl(iters = 120L, print = 0L))))
    expect_s3_class(fit, "nlmixr2FitData")
    expect_true(is.finite(fit$objf))
    expect_true(all(c("IPRED", "CWRES") %in% names(fit)))
    ## population parameter table present with finite estimates
    expect_true(is.data.frame(fit$parFixedDf))
    expect_true(all(is.finite(fit$parFixedDf$Estimate)))
    ## ADVI artifacts carried on the fit env
    expect_false(is.null(fit$env$adviElbo))
    expect_false(is.null(fit$env$adviState))
    ## the optimization walk is standard parHistData (captured even with print=0)
    .ph <- fit$parHist
    expect_true(is.data.frame(.ph))
    expect_true(all(c("iter", "tka", "add.sd", "o(eta.ka)") %in% names(.ph)))
    expect_gte(max(.ph$iter), 120L)
    ## the finalize reuses the loop's compiled models (no symengine rebuild)
    expect_false(is.null(fit$env$foceiModel))
  })

  test_that("adviControl(tol=) stops early on ELBO convergence", {
    ## The mechanism must be OBSERVABLE, not just "the fit still works": with a
    ## loose tolerance the loop must stop before `iters`, and with tol=0 it must
    ## run every iteration.  Before this was wired up, `tol` was documented but
    ## never reached the C++ loop, so both runs were identical -- a test that
    ## only checked the fit succeeded would not have caught that.
    .base <- function(tol) adviControl(iters = 400L, seed = 7L, print = 0L,
                                       returnAdvi = TRUE, tol = tol,
                                       evalElbo = 25L)
    rOff <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "advi", control = .base(0))))
    expect_equal(rOff$itRun, 400L)
    expect_false(isTRUE(rOff$tolStopped))
    expect_equal(length(rOff$elbo), 400L)

    rOn <- suppressMessages(suppressWarnings(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "advi", control = .base(1))))
    expect_true(isTRUE(rOn$tolStopped))      # rel change < 1 always -> first check
    expect_lt(rOn$itRun, 400L)
    ## the reported trace is truncated to what actually ran, not zero-padded
    expect_equal(length(rOn$elbo), rOn$itRun)
    expect_true(all(is.finite(rOn$elbo)))
  })
})
