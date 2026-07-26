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
    f <- nlmixr2(theo, nlmixr2data::theo_sd, est = "vae",
                 control = vaeControl(itersBurnIn = 80L, klWarmup = 40L, gammaIter = 120L,
                                      iters = 160L, hiddenDim = 25L, seed = 1L))
    expect_s3_class(f, "nlmixr2FitData")

    ## selected covariate coefficients are population parameters in the fit.
    ## A continuous coefficient is named beta_<par>_<cov>_<shape> and which shape
    ## family wins is data-driven, so match the covariate and read the shape back.
    .pf <- if (is.null(rownames(f$parFixed))) f$parFixed$Parameter else rownames(f$parFixed)
    .beta <- function(par) grep(paste0("^beta_", par, "_WT(_|$)"), .pf, value = TRUE)
    expect_length(.beta("lka"), 1L)
    expect_length(.beta("lV"), 1L)
    ## objDf / objective present, tagged est = vae
    expect_true(is.finite(f$objDf$OBJF[1]))

    ## FINAL model carries the exact centered covariate expression -- the one
    ## belonging to the shape each coefficient was named for
    .fin <- vapply(f$finalUi$lstExpr, function(e) deparse1(e), character(1))
    .term <- function(nm) {
      switch(sub("^beta_[^_]+_WT_?", "", nm),
             power = "log\\(WT/[0-9.]+\\)",
             log = "log\\(WT\\)",
             lin = "\\(WT - [0-9.]+\\)",
             identity = "WT",
             center = "\\(WT/[0-9.]+\\)",
             stop("unexpected coefficient name: ", nm))
    }
    for (.p in c("lka", "lV")) {
      .b <- .beta(.p)
      expect_true(any(grepl(paste0(.b, " \\* ", .term(.b)), .fin)))
    }
    expect_true(any(grepl("^ke <- exp\\(lke \\+ eta.ke\\)$", .fin)))

    ## ORIGINAL model recoverable via $uiIni (structure differs -- no covariates)
    .ini <- vapply(f$uiIni$lstExpr, function(e) deparse1(e), character(1))
    expect_false(any(grepl("beta_lka_WT", .ini)))
    expect_true(any(grepl("^ka <- exp\\(lka \\+ eta.ka\\)$", .ini)))
    ## $iniDf0 is the ORIGINAL model's iniDf (no covariate thetas)
    expect_false(any(grepl("^beta_", f$iniDf0$name)))
    expect_true(all(c("lka", "lke", "lV") %in% f$iniDf0$name))

    ## EBEs from the encoder
    expect_true(all(c("eta.ka", "eta.ke", "eta.V") %in% names(f$eta)))
  })
})
