## Cross-check est="advi" against the neonatal weight-progression case study
## (Rohleff et al. 2025) used for est="vae": a 5-parameter turnover growth model
## (log-linked W0/kin/TL/koutmax/T50, an ODE, W(0)=W0, combined error) over 189
## neonates.  ADVI is run WITHOUT covariate selection (unlike the vae example),
## so it fits the base structural model; the population birth weight must land in
## the same physiological range the vae fit reports.  Weekly batch; skipped when
## the dataset is unavailable in the installed nlmixr2data.

nmTest({
  test_that("est='advi' fits the neonatal turnover model (no covariate selection)", {
    skip_on_cran()
    skip_if_not(exists("neonatal_wt", where = asNamespace("nlmixr2data")),
                "nlmixr2data::neonatal_wt not available")
    neonatal <- function() {
      ini({ lW0 <- log(3000); lkin <- log(30); lTL <- log(2)
        lkoutmax <- log(0.05); lT50 <- log(1)
        eta.W0 ~ 0.001; eta.kin ~ 0.01; eta.TL ~ 0.05
        eta.koutmax ~ 0.05; eta.T50 ~ 0.05
        add.err <- 30; prop.err <- 0.02 })
      model({ W0 <- exp(lW0 + eta.W0); kin <- exp(lkin + eta.kin); TL <- exp(lTL + eta.TL)
        koutmax <- exp(lkoutmax + eta.koutmax); T50 <- exp(lT50 + eta.T50)
        kprod <- kin * expit(2 * (t - TL)); kelim <- koutmax * (1 - t / (T50 + t))
        d/dt(W) <- kprod - kelim * W; W(0) <- W0
        W ~ add(add.err) + prop(prop.err) })
    }
    dat <- get("neonatal_wt", envir = asNamespace("nlmixr2data"))
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(neonatal, dat, est = "advi",
              control = adviControl(iters = 400L, print = 0L, returnAdvi = TRUE))))

    ## all five structural typical values + the two error params estimated finite
    expect_true(all(is.finite(fit$theta)))
    expect_true(all(is.finite(fit$popOmega)) && all(fit$popOmega > 0))
    ## typical birth weight in a physiologically sane range (grams), matching the
    ## vae fit's canonical result
    bw <- exp(fit$theta[1])
    expect_gt(bw, 2500); expect_lt(bw, 5000)
  })

  test_that("est='advi' assembles a full neonatal fit object", {
    skip_on_cran()
    skip_if_not(exists("neonatal_wt", where = asNamespace("nlmixr2data")),
                "nlmixr2data::neonatal_wt not available")
    neonatal <- function() {
      ini({ lW0 <- log(3000); lkin <- log(30); lTL <- log(2)
        lkoutmax <- log(0.05); lT50 <- log(1)
        eta.W0 ~ 0.001; eta.kin ~ 0.01; eta.TL ~ 0.05
        eta.koutmax ~ 0.05; eta.T50 ~ 0.05
        add.err <- 30; prop.err <- 0.02 })
      model({ W0 <- exp(lW0 + eta.W0); kin <- exp(lkin + eta.kin); TL <- exp(lTL + eta.TL)
        koutmax <- exp(lkoutmax + eta.koutmax); T50 <- exp(lT50 + eta.T50)
        kprod <- kin * expit(2 * (t - TL)); kelim <- koutmax * (1 - t / (T50 + t))
        d/dt(W) <- kprod - kelim * W; W(0) <- W0
        W ~ add(add.err) + prop(prop.err) })
    }
    dat <- get("neonatal_wt", envir = asNamespace("nlmixr2data"))
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(neonatal, dat, est = "advi", control = adviControl(iters = 400L, print = 0L))))
    expect_s3_class(fit, "nlmixr2FitData")
    expect_true(is.finite(fit$objf))
    expect_true(all(c("IPRED", "CWRES") %in% names(fit)))
  })
})
