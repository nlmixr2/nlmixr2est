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

  test_that("nonMuTheta='fix' holds the injected omega fixed", {
    ui <- rxode2::assertRxUi(.nmt())
    r <- suppressWarnings(suppressMessages(
      .preProcessVaeNonMuTheta(ui, "vae", NULL, vaeControl(nonMuTheta = "fix", nonMuEtaOmega = 0.001))))
    idf <- r$ui$iniDf
    inj <- idf[!is.na(idf$neta1) & idf$neta1 == idf$neta2 & idf$name != "eta.kout", ]
    expect_true(all(inj$est == 0.001))
    expect_true(all(inj$fix))
  })

  test_that("nonMuTheta='none' injects nothing", {
    ui <- rxode2::assertRxUi(.nmt())
    expect_null(.preProcessVaeNonMuTheta(ui, "vae", NULL, vaeControl(nonMuTheta = "none")))
  })

  test_that("vaeControl exposes nonMuTheta/nonMuEtaOmega with the documented defaults", {
    expect_equal(vaeControl()$nonMuTheta, "eta")
    expect_equal(vaeControl()$nonMuEtaOmega, 0.01)
    expect_equal(vaeControl(nonMuTheta = "fix")$nonMuTheta, "fix")
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
                      nGradStep = 2L, covariateSelection = FALSE, print = 0L)
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
})
