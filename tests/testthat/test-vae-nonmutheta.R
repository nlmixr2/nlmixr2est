## est="vae" only estimates parameters in the latent space (thetas with a random
## effect).  nonMuTheta injects a small eta into each non-mu-referenced structural
## theta so it is estimated, reports it as theta+mean(eta), and drops the temporary
## eta.  Unit tests (detection/injection/control) are fast; one small end-to-end fit
## confirms the parameter actually moves and the eta is dropped.

nmTest({
  .nmt <- function() {
    ini({ tR0 <- log(100); tkout <- log(1); tIC50 <- log(5); eta.kout ~ 1; prop.err <- 0.3 })
    model({ R0 <- exp(tR0); kout <- exp(tkout + eta.kout); IC50 <- exp(tIC50)
      kel <- 0.1; V <- 10; Imax <- 0.9; kin <- R0 * kout; R(0) <- R0
      d/dt(A) <- -kel * A; d/dt(R) <- kin * (1 - Imax * (A / V) / (IC50 + (A / V))) - kout * R
      R ~ prop(prop.err) })
  }

  test_that("non-mu-referenced structural thetas are detected", {
    ui <- rxode2::assertRxUi(.nmt())
    expect_setequal(.vaeNonMuThetas(ui), c("tR0", "tIC50"))
  })

  test_that("nonMuTheta='eta' injects estimated etas and mu-references the thetas", {
    ui <- rxode2::assertRxUi(.nmt())
    r <- suppressWarnings(suppressMessages(
      .preProcessVaeNonMuTheta(ui, "vae", NULL, vaeControl(nonMuTheta = "eta", nonMuEtaOmega = 0.02))))
    idf <- r$ui$iniDf
    inj <- idf[!is.na(idf$neta1) & idf$neta1 == idf$neta2 & idf$name != "eta.kout", ]
    expect_equal(nrow(inj), 2L)                 # eta for tR0 and tIC50
    expect_true(all(inj$est == 0.02))           # nonMuEtaOmega
    expect_true(all(!inj$fix))                  # estimated (not fixed)
    expect_true(all(c("tR0", "tIC50") %in% r$ui$muRefDataFrame$theta))  # now mu-referenced
  })

  test_that("nonMuTheta='fix' holds the injected omega AND the theta fixed", {
    ui <- rxode2::assertRxUi(.nmt())
    r <- suppressWarnings(suppressMessages(
      .preProcessVaeNonMuTheta(ui, "vae", NULL, vaeControl(nonMuTheta = "fix", nonMuEtaOmega = 0.001))))
    idf <- r$ui$iniDf
    inj <- idf[!is.na(idf$neta1) & idf$neta1 == idf$neta2 & idf$name != "eta.kout", ]
    expect_true(all(inj$est == 0.001))
    expect_true(all(inj$fix))
    ## the paired structural thetas are held fixed too (the fixed effect is not
    ## estimated in "fix" mode)
    expect_true(all(idf$fix[idf$name %in% c("tR0", "tIC50") & !is.na(idf$ntheta)]))
    ## prep flags those latent dims as zPopFix so the M-step holds them at ini
    p <- .vaeDataPrep(r$ui, data.frame(ID = 1L, TIME = c(0, 1), DV = c(1, 2), AMT = c(1, 0)),
                      vaeControl(nonMuTheta = "fix"))
    expect_true(any(p$zPopFix))
  })

  test_that("nonMuTheta='none' injects nothing", {
    ui <- rxode2::assertRxUi(.nmt())
    expect_null(.preProcessVaeNonMuTheta(ui, "vae", NULL, vaeControl(nonMuTheta = "none")))
  })

  test_that("nonMuTheta='regress' injects no eta and leaves the fixed thetas intact", {
    ui <- rxode2::assertRxUi(.nmt())
    ## the hook makes no model change (the thetas are regressed later in the M-step)
    r <- suppressWarnings(suppressMessages(
      .preProcessVaeNonMuTheta(ui, "vae", NULL, vaeControl(nonMuTheta = "regress"))))
    expect_null(r)
    ## no eta injected: the UI still carries only the original eta.kout
    expect_equal(ui$iniDf[!is.na(ui$iniDf$neta1) & ui$iniDf$neta1 == ui$iniDf$neta2, "name"],
                 "eta.kout")
    ## prep exposes the regress-target theta indices + ini bounds
    p <- .vaeDataPrep(ui, data.frame(ID = 1L, TIME = c(0, 1), DV = c(1, 2), AMT = c(1, 0)),
                      vaeControl(nonMuTheta = "regress"))
    expect_setequal(p$regressNames, c("tR0", "tIC50"))
    expect_equal(length(p$regressThetaIdx0), 2L)
    ## "eta" mode leaves the regress fields empty
    p2 <- .vaeDataPrep(ui, data.frame(ID = 1L, TIME = c(0, 1), DV = c(1, 2), AMT = c(1, 0)),
                       vaeControl(nonMuTheta = "eta"))
    expect_equal(length(p2$regressNames), 0L)
  })

  test_that("vaeControl exposes nonMuTheta/nonMuEtaOmega with the documented defaults", {
    expect_equal(vaeControl()$nonMuTheta, "regress")
    expect_equal(vaeControl()$nonMuEtaOmega, 0.01)
    expect_equal(vaeControl(nonMuTheta = "fix")$nonMuTheta, "fix")
    expect_equal(vaeControl(nonMuTheta = "eta")$nonMuTheta, "eta")
    expect_error(vaeControl(nonMuTheta = "bogus"))
    expect_error(vaeControl(nonMuEtaOmega = -1))
  })

  test_that("est=vae estimates a non-mu theta, drops the eta, and warns", {
    skip_on_cran()
    theo <- function() {
      ini({ lka <- log(1.8); lke <- log(0.086); lV <- log(32); eta.ka ~ 0.3; add.err <- 0.7 })
      model({ ka <- exp(lka + eta.ka); ke <- exp(lke); V <- exp(lV)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ctl <- vaeControl(itersBurnIn = 5L, iters = 12L, klWarmup = 5L, gammaIter = 8L,
                      nGradStep = 2L, covariateSelection = FALSE, print = 0L,
                      nonMuTheta = "eta")
    fit <- suppressWarnings(suppressMessages(
      nlmixr2(theo, nlmixr2data::theo_sd, est = "vae", control = ctl)))
    ## lke and lV (no eta in the input) are now estimated (moved off their ini)
    expect_gt(abs(fit$theta[["lke"]] - log(0.086)), 1e-4)
    expect_gt(abs(fit$theta[["lV"]] - log(32)), 1e-4)
    ## the injected etas are dropped from the reported model (only eta.ka remains)
    expect_equal(fit$iniDf[!is.na(fit$iniDf$neta1) & fit$iniDf$neta1 == fit$iniDf$neta2, "name"],
                 "eta.ka")
    ## the warning is recorded in $runInfo
    expect_true(any(grepl("non-mu-referenced", fit$runInfo)))

    ## nonMuTheta="none" keeps them frozen at their ini values
    fit0 <- suppressWarnings(suppressMessages(
      nlmixr2(theo, nlmixr2data::theo_sd, est = "vae",
              control = vaeControl(itersBurnIn = 5L, iters = 12L, klWarmup = 5L, gammaIter = 8L,
                                   nGradStep = 2L, covariateSelection = FALSE, print = 0L,
                                   nonMuTheta = "none"))))
    expect_equal(fit0$theta[["lke"]], log(0.086), tolerance = 1e-6)
    expect_equal(fit0$theta[["lV"]], log(32), tolerance = 1e-6)
  })

  test_that("nonMuTheta='regress' estimates a non-mu theta by bounded bobyqa, no eta", {
    skip_on_cran()
    theo <- function() {
      ini({ lka <- log(1.8); lke <- c(log(0.01), log(0.086), log(1)); lV <- log(32)
        eta.ka ~ 0.3; add.err <- 0.7 })
      model({ ka <- exp(lka + eta.ka); ke <- exp(lke); V <- exp(lV)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ctl <- vaeControl(itersBurnIn = 5L, iters = 15L, klWarmup = 5L, gammaIter = 10L,
                      nGradStep = 2L, covariateSelection = FALSE, print = 0L,
                      nonMuTheta = "regress")
    fit <- suppressWarnings(suppressMessages(
      nlmixr2(theo, nlmixr2data::theo_sd, est = "vae", control = ctl)))
    ## the regression moved the non-mu thetas off their ini (evidence the M-step
    ## bobyqa step ran) ...
    expect_gt(abs(fit$theta[["lke"]] - log(0.086)), 1e-4)
    expect_gt(abs(fit$theta[["lV"]] - log(32)), 1e-4)
    ## ... and stayed within the ini() bounds for the bounded parameter
    expect_gte(fit$theta[["lke"]], log(0.01))
    expect_lte(fit$theta[["lke"]], log(1))
    ## NO eta was injected: only the real eta.ka remains
    expect_equal(fit$iniDf[!is.na(fit$iniDf$neta1) & fit$iniDf$neta1 == fit$iniDf$neta2, "name"],
                 "eta.ka")
    ## the regress note is recorded in $runInfo
    expect_true(any(grepl("bobyqa regression", fit$runInfo)))
    ## the regressed thetas are surfaced in the iteration-print / parameter history
    ## (not just estimated silently) and their final walk value matches the fit
    expect_true(all(c("lke", "lV") %in% names(fit$parHistData)))
    .last <- fit$parHistData[nrow(fit$parHistData), ]
    expect_equal(.last[["lke"]], exp(fit$theta[["lke"]]), tolerance = 1e-3)
    expect_equal(.last[["lV"]], exp(fit$theta[["lV"]]), tolerance = 1e-3)
  })

  test_that("nonMuTheta='fix' holds the theta at ini and omits it from the print", {
    skip_on_cran()
    theo <- function() {
      ini({ lka <- log(1.8); lke <- log(0.086); lV <- log(32); eta.ka ~ 0.3; add.err <- 0.7 })
      model({ ka <- exp(lka + eta.ka); ke <- exp(lke); V <- exp(lV)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ctl <- vaeControl(itersBurnIn = 5L, iters = 12L, klWarmup = 5L, gammaIter = 8L,
                      nGradStep = 2L, covariateSelection = FALSE, print = 0L,
                      nonMuTheta = "fix")
    fit <- suppressWarnings(suppressMessages(
      nlmixr2(theo, nlmixr2data::theo_sd, est = "vae", control = ctl)))
    ## the fixed effects stay exactly at their ini() value (not estimated) ...
    expect_equal(fit$theta[["lke"]], log(0.086), tolerance = 1e-8)
    expect_equal(fit$theta[["lV"]], log(32), tolerance = 1e-8)
    ## ... are reported as fixed, with the injected eta dropped ...
    idf <- fit$iniDf
    expect_true(all(idf$fix[idf$name %in% c("lke", "lV") & !is.na(idf$ntheta)]))
    expect_equal(idf[!is.na(idf$neta1) & idf$neta1 == idf$neta2, "name"], "eta.ka")
    ## ... and are NOT shown in the iteration table / parameter history
    expect_false(any(c("lke", "lV") %in% names(fit$parHistData)))
  })
})
