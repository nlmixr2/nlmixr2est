# nlme covMethod: the default recomputes the FOCEI observed-information
# covariance at the converged nlme estimates, keeping nlme's own covariance
# selectable via the "nlme" token.  Weekly batch (multi-iteration fits).

nmTest({
  mod <- function() {
    ini({
      tka <- 0.45; tcl <- log(c(0, 2.7, 100)); tv <- 3.45
      eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  test_that("nlme default covMethod recomputes the analytic covariance", {
    skip_on_cran()
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(mod, nlmixr2data::theo_sd, est = "nlme", control = nlmeControl(print = 0))))
    ## analytic (or its FD fallback), never the legacy nlme label by default
    expect_false(identical(fit$covMethod, "nlme"))
    expect_true(all(is.finite(fit$parFixedDf$SE)))
    expect_true(all(fit$parFixedDf$SE > 0))
    ## nlme's own covariance is retained and selectable
    expect_true("nlme" %in% names(fit$env$covList))
    expect_error(setCov(fit, "nlme"), NA)
    expect_identical(fit$covMethod, "nlme")
  })

  test_that("covMethod='nlme' keeps nlme's own covariance (no recompute)", {
    skip_on_cran()
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(mod, nlmixr2data::theo_sd, est = "nlme",
              control = nlmeControl(print = 0, covMethod = "nlme"))))
    expect_identical(fit$covMethod, "nlme")
  })
})
