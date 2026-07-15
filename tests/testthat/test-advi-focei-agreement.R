## Point-estimate ADVI recovers the maximum-likelihood estimates: on a model with
## mu-referenced structural thetas (tka, tcl), a non-mu structural theta (tv,
## estimated through the theta-sensitivity gradient), and a residual sigma, the
## ADVI population estimates agree with FOCEI within tolerance.  Weekly batch
## (multi-iteration fit).

nmTest({
  test_that("point-estimate ADVI agrees with FOCEI on theo_sd", {
    skip_on_cran()
    mod <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; eta.cl ~ 0.3; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka*depot; d/dt(center) <- ka*depot - cl/v*center
        cp <- center/v; cp ~ add(add.sd) })
    }
    fF <- suppressMessages(suppressWarnings(
      nlmixr2(mod, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(print = 0L))))
    fA <- suppressMessages(suppressWarnings(
      nlmixr2(mod, nlmixr2data::theo_sd, est = "advi",
              control = adviControl(iters = 600L, print = 0L, returnAdvi = TRUE))))

    ## typical values (log scale): mu-ref tka/tcl, non-mu tv, sigma add.sd
    expect_equal(unname(fA$theta), unname(fF$theta), tolerance = 0.1)
    ## between-subject variances agree to ~30% (mean-field is approximate)
    expect_equal(unname(fA$popOmega), unname(diag(fF$omega)), tolerance = 0.3)
  })

  test_that("point-estimate ADVI builds a fit whose objf is near FOCEI's", {
    skip_on_cran()
    mod <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; eta.cl ~ 0.3; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka*depot; d/dt(center) <- ka*depot - cl/v*center
        cp <- center/v; cp ~ add(add.sd) })
    }
    fF <- suppressMessages(suppressWarnings(
      nlmixr2(mod, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(print = 0L))))
    fA <- suppressMessages(suppressWarnings(
      nlmixr2(mod, nlmixr2data::theo_sd, est = "advi",
              control = adviControl(iters = 600L, print = 0L))))
    expect_s3_class(fA, "nlmixr2FitData")
    ## the FOCEi objective at the ADVI estimates is within a few points of FOCEI's
    expect_lt(abs(fA$objf - fF$objf), 15)
  })
})
