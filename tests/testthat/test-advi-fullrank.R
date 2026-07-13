## Block full-rank ADVI (adviFamily="fullRank"): each subject gets a full
## neta x neta Cholesky factor L_i, so q(eta_i)=N(mu_i, L_i L_i^T) captures
## within-subject posterior correlation the mean-field family cannot.  It
## recovers the FOCEi ML population estimates, and its per-subject posterior
## carries non-trivial off-diagonal covariance.  Weekly batch (multi-iteration).

nmTest({
  mod <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; eta.cl ~ 0.3; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      d/dt(depot) <- -ka*depot; d/dt(center) <- ka*depot - cl/v*center
      cp <- center/v; cp ~ add(add.sd) })
  }

  test_that("full-rank ADVI recovers the FOCEI estimates", {
    skip_on_cran()
    fF <- suppressMessages(suppressWarnings(
      nlmixr2(mod, nlmixr2data::theo_sd, est = "focei", control = foceiControl(print = 0L))))
    fA <- suppressMessages(suppressWarnings(
      nlmixr2(mod, nlmixr2data::theo_sd, est = "advi",
              control = adviControl(iters = 500L, print = 0L, returnAdvi = TRUE,
                                    adviFamily = "fullRank"))))
    expect_identical(fA$family, "fullRank")
    expect_equal(unname(fA$theta), unname(fF$theta), tolerance = 0.1)
    expect_equal(unname(fA$popOmega), unname(diag(fF$omega)), tolerance = 0.3)
  })

  test_that("full-rank ADVI captures within-subject posterior correlation", {
    skip_on_cran()
    fA <- suppressMessages(suppressWarnings(
      nlmixr2(mod, nlmixr2data::theo_sd, est = "advi",
              control = adviControl(iters = 500L, print = 0L, returnAdvi = TRUE,
                                    adviFamily = "fullRank"))))
    ## per-subject L is packed lower-tri (neta=2 -> 3 cols: L11, L21, L22).  A
    ## non-negligible off-diagonal (L21) for at least some subjects means the
    ## posterior covariance L L^T has real off-diagonal structure.
    Lp <- fA$scale
    expect_equal(ncol(Lp), 3L)
    expect_gt(max(abs(Lp[, 2L])), 1e-3)
  })

  test_that("full-rank ADVI assembles a full nlmixr2FitData", {
    skip_on_cran()
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(mod, nlmixr2data::theo_sd, est = "advi",
              control = adviControl(iters = 400L, print = 0L, adviFamily = "fullRank"))))
    expect_s3_class(fit, "nlmixr2FitData")
    expect_true(is.finite(fit$objf))
    expect_true(all(c("IPRED", "CWRES") %in% names(fit)))
  })
})
