## Phase 6: est="vae" returns a standard nlmixr2FitData assembled WITHOUT running
## focei estimation. The final model carries the selected covariate effects
## (exact centered expressions); the original model is recoverable via $uiIni /
## $iniDf0 (the structure changes under covariate selection); the objective is
## the VAE linearization -2LL; the EBEs come from the encoder. Slow; CRAN-skipped.

nmTest({
  test_that("est=\"vae\" builds a nlmixr2FitData with original/final model accessors", {
    skip_on_cran()
    theo <- function() {
      ini({
        lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03
        add.err <- 0.7
      })
      model({
        ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot
        d/dt(central) = ka * depot - ke * central
        cp <- central / V
        cp ~ add(add.err)
      })
    }
    ## shapes/covCenterType pinned as in test-vae-covariate-selection.R: the
    ## expanded shape
    ## search would otherwise decide the coefficient NAMES (lka picks "lin" and
    ## lV "power" by default), making this test about which shape wins rather
    ## than about the fit accessors it is here to check.
    f <- nlmixr2(theo, nlmixr2data::theo_sd, est = "vae",
                 control = vaeControl(itersBurnIn = 80L, klWarmup = 40L, gammaIter = 120L,
                                      iters = 160L, hiddenDim = 25L, seed = 1L,
                                      shapes = "power", covCenterType = "mean"))
    expect_s3_class(f, "nlmixr2FitData")

    ## selected covariate coefficients are population parameters in the fit
    .pf <- if (is.null(rownames(f$parFixed))) f$parFixed$Parameter else rownames(f$parFixed)
    expect_true(all(c("beta.lka.WT.power", "beta.lV.WT.power") %in% .pf))
    ## objDf / objective present, tagged est = vae
    expect_true(is.finite(f$objDf$OBJF[1]))

    ## FINAL model carries the exact centered covariate expression
    .fin <- vapply(f$finalUi$lstExpr, function(e) deparse1(e), character(1))
    expect_true(any(grepl("beta.lka.WT.power \\* log\\(WT/", .fin, fixed = FALSE)))
    expect_true(any(grepl("beta.lV.WT.power \\* log\\(WT/", .fin, fixed = FALSE)))
    expect_true(any(grepl("^ke <- exp\\(lke \\+ eta.ke\\)$", .fin)))

    ## ORIGINAL model recoverable via $uiIni (structure differs -- no covariates)
    .ini <- vapply(f$uiIni$lstExpr, function(e) deparse1(e), character(1))
    expect_false(any(grepl("beta.lka.WT", .ini, fixed = TRUE)))
    expect_true(any(grepl("^ka <- exp\\(lka \\+ eta.ka\\)$", .ini)))
    ## $iniDf0 is the ORIGINAL model's iniDf (no covariate thetas)
    expect_false(any(grepl("^beta\\.", f$iniDf0$name)))
    expect_true(all(c("lka", "lke", "lV") %in% f$iniDf0$name))

    ## EBEs from the encoder
    expect_true(all(c("eta.ka", "eta.ke", "eta.V") %in% names(f$eta)))
  })
})
