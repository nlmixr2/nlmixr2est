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
  })
})
