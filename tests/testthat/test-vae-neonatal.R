## Phase 8 (Milestone C): the neonatal weight-progression case study from Rohleff
## et al. (2025). A turnover growth model (5 log-linked params W0/kin/TL/koutmax/
## T50, an ODE, W(0)=W0, combined error) with 5 candidate covariates, fit by the
## inner-driver VAE with BICc-ELBO covariate selection over 189 neonates. The
## canonical result is that birth weight (W0) depends on gestational age (GA).

nmTest({
  test_that("est=vae fits the neonatal weight model and selects GA on birth weight", {
    skip_on_cran()
    neonatal <- function() {
      ini({
        lW0 <- log(3000); lkin <- log(30); lTL <- log(2); lkoutmax <- log(0.05); lT50 <- log(1)
        eta.W0 ~ 0.001; eta.kin ~ 0.01; eta.TL ~ 0.05; eta.koutmax ~ 0.05; eta.T50 ~ 0.05
        add.err <- 30; prop.err <- 0.02
      })
      model({
        W0 <- exp(lW0 + eta.W0); kin <- exp(lkin + eta.kin); TL <- exp(lTL + eta.TL)
        koutmax <- exp(lkoutmax + eta.koutmax); T50 <- exp(lT50 + eta.T50)
        kprod <- kin * expit(2 * (t - TL))
        kelim <- koutmax * (1 - t / (T50 + t))
        d/dt(W) <- kprod - kelim * W
        W(0) <- W0
        W ~ add(add.err) + prop(prop.err)
      })
    }
    dat <- nlmixr2data::neonatal_wt
    ctl <- vaeControl(itersBurnIn = 60L, iters = 120L, klWarmup = 40L, gammaIter = 90L,
                      nGradStep = 4L, covariateSelection = TRUE, print = 0L,
                      sigma0 = c(1e-3, 1e-2, 1e-1, 1e-1, 1e-1), returnVae = TRUE)
    fit <- suppressMessages(suppressWarnings(nlmixr2(neonatal, dat, est = "vae", control = ctl)))

    ## all five structural params + the (combined) error estimated finite
    expect_length(fit$zPop, 5L)
    expect_true(all(is.finite(fit$zPop)))
    expect_true(all(is.finite(fit$a)) && all(fit$a >= 0))
    ## typical birth weight in a physiologically sane range (grams)
    expect_gt(exp(fit$zPop[1]), 2500)
    expect_lt(exp(fit$zPop[1]), 5000)

    ## covariate selection: gestational age (GA) drives birth weight (W0) -- the
    ## strongest, canonical effect (Fig 4, VAE column)
    dimnames(fit$selected) <- list(c("W0", "kin", "TL", "koutmax", "T50"), fit$covNames)
    expect_true(fit$selected["W0", "GA"])

    ## The reference M-step objective must reproduce this fit EXACTLY.  Every
    ## structural parameter here is mu-referenced, so the model has no non-mu
    ## theta for mStepObjective= to act on -- the option is out of the picture and
    ## the two runs have to agree bit for bit.  That pins two things at once: the
    ## deviation is confined to the non-mu theta M-step (it does NOT rescore the
    ## covariate search or the encoder ELBO), and the covariate set this case
    ## study selects is therefore NOT attributable to the objective choice.
    ctlRef <- ctl
    ctlRef$mStepObjective <- "elbo"
    fitRef <- suppressMessages(suppressWarnings(
      nlmixr2(neonatal, dat, est = "vae", control = ctlRef)))
    expect_identical(fitRef$zPop, fit$zPop)
    expect_identical(fitRef$selected, fit$selected)
    expect_identical(fitRef$omega, fit$omega)
    expect_identical(fitRef$a, fit$a)
  })

  test_that("est=vae builds a full neonatal fit object (objective + tables)", {
    skip_on_cran()
    skip_if_not(identical(Sys.getenv("NLMIXR2_VAE_SLOW"), "true"),
                "set NLMIXR2_VAE_SLOW=true to run the full neonatal fit assembly")
    neonatal <- function() {
      ini({
        lW0 <- log(3000); lkin <- log(30); lTL <- log(2); lkoutmax <- log(0.05); lT50 <- log(1)
        eta.W0 ~ 0.001; eta.kin ~ 0.01; eta.TL ~ 0.05; eta.koutmax ~ 0.05; eta.T50 ~ 0.05
        add.err <- 30; prop.err <- 0.02
      })
      model({
        W0 <- exp(lW0 + eta.W0); kin <- exp(lkin + eta.kin); TL <- exp(lTL + eta.TL)
        koutmax <- exp(lkoutmax + eta.koutmax); T50 <- exp(lT50 + eta.T50)
        kprod <- kin * expit(2 * (t - TL)); kelim <- koutmax * (1 - t / (T50 + t))
        d/dt(W) <- kprod - kelim * W; W(0) <- W0
        W ~ add(add.err) + prop(prop.err)
      })
    }
    ctl <- vaeControl(itersBurnIn = 60L, iters = 120L, klWarmup = 40L, gammaIter = 90L,
                      nGradStep = 4L, covariateSelection = TRUE, print = 0L,
                      sigma0 = c(1e-3, 1e-2, 1e-1, 1e-1, 1e-1))
    fit <- suppressMessages(suppressWarnings(nlmixr2(neonatal, nlmixr2data::neonatal_wt, est = "vae", control = ctl)))
    expect_s3_class(fit, "nlmixr2FitData")
    expect_true(is.finite(fit$objf))
    expect_true(all(c("CWRES", "IPRED") %in% names(fit)))
  })
})
