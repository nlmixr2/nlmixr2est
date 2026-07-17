## Full-Bayes ADVI (pointEstimate=FALSE): the free population parameters get a
## full-rank variational Gaussian q(phi)=N(mPop, Lpop Lpop^T) on top of the
## per-subject q(eta_i).  The population posterior MEANS recover the FOCEi ML
## estimates, and the population variational covariance (Lpop Lpop^T) gives
## sensible parameter uncertainties (not the inflated flat-direction variance).
## Works for both per-subject families.  Weekly batch (multi-iteration).

nmTest({
  mod <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; eta.cl ~ 0.3; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      d/dt(depot) <- -ka*depot; d/dt(center) <- ka*depot - cl/v*center
      cp <- center/v; cp ~ add(add.sd) })
  }

  for (fam in c("meanField", "fullRank")) {
    test_that(paste0("full-Bayes ADVI (", fam, ") recovers FOCEI with sane uncertainties"), {
      skip_on_cran()
      fF <- suppressMessages(suppressWarnings(
        nlmixr2(mod, nlmixr2data::theo_sd, est = "focei", control = foceiControl(print = 0L))))
      fA <- suppressMessages(suppressWarnings(
        nlmixr2(mod, nlmixr2data::theo_sd, est = "advi",
                control = adviControl(iters = 500L, print = 0L, returnAdvi = TRUE,
                                      pointEstimate = FALSE, adviFamily = fam))))
      expect_false(fA$pointEstimate)
      ## posterior MEANS agree with FOCEI
      expect_equal(unname(fA$theta), unname(fF$theta), tolerance = 0.15)
      expect_equal(unname(fA$popOmega), unname(diag(fF$omega)), tolerance = 0.3)
      ## population variational covariance: symmetric PSD, sane theta SEs
      cov <- fA$adviCov
      expect_equal(dim(cov), rep(length(fA$prep$theta) + fA$prep$neta -
                                 sum(fA$prep$thetaFix) - sum(fA$prep$omegaFix), 2L)[c(1, 1)])
      sds <- sqrt(diag(cov))
      expect_true(all(is.finite(sds)) && all(sds > 0))
      ## the typical-value SEs are small (a few percent), NOT the inflated
      ## flat-direction variance (~2) the mu-referencing would give if unhandled
      expect_lt(max(sds[seq_len(sum(!fA$prep$thetaFix))]), 0.5)
    })
  }

  test_that("full-Bayes ADVI assembles a full nlmixr2FitData", {
    skip_on_cran()
    ## the default covMethod is "advi": the population variational covariance is
    ## the fit's SE source
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(mod, nlmixr2data::theo_sd, est = "advi",
              control = adviControl(iters = 400L, print = 0L, pointEstimate = FALSE))))
    expect_s3_class(fit, "nlmixr2FitData")
    expect_true(is.finite(fit$objf))
    expect_true(all(c("IPRED", "CWRES") %in% names(fit)))
    expect_false(is.null(fit$env$adviCov))
    ## the population variational covariance is the fit's SE source
    expect_identical(fit$covMethod, "advi")
    expect_true(all(is.finite(fit$parFixedDf$SE)))
    expect_true(all(fit$parFixedDf$SE > 0))
  })

  test_that("full-Bayes ADVI covMethod='analytic' is not clobbered to advi", {
    skip_on_cran()
    ## an explicit non-"advi" covMethod runs the FOCEi covariance chain on the
    ## full inner model; the full-Bayes path must NOT overwrite it with the
    ## variational covariance
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(mod, nlmixr2data::theo_sd, est = "advi",
              control = adviControl(iters = 400L, print = 0L, pointEstimate = FALSE,
                                    covMethod = "analytic"))))
    expect_false(identical(fit$covMethod, "advi"))
    expect_false(is.null(fit$env$adviCov))          # the variational cov is still kept as an artifact
    expect_true(all(is.finite(fit$parFixedDf$SE)))
    expect_true(all(fit$parFixedDf$SE > 0))
  })
})
