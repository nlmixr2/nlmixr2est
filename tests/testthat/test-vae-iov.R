## IOV support for est="vae". Enabled purely by the "iov" attribute on
## nlmixr2Est.vae: the .uiApplyIov hook materializes each occasion-level random
## effect (iov.x ~ v | occ) into per-occasion fixed-variance(1) etas, which the
## VAE handles as non-mu-referenced (theta=0-centered) etas with held omega.
## theta+eta+iov (lka+eta.ka+iov.ka) is already characterized as a mu-referenced
## theta+eta expression by rxode2 (the id-level eta pairs with the theta), so no
## rxode2 change is needed.

nmTest({
  .vaeIovMod <- function() {
    ini({ lka <- log(1.8); lke <- log(0.086); lV <- log(32)
      eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03
      iov.ka ~ 0.1 | occ
      add.err <- 0.7 })
    model({ ka <- exp(lka + eta.ka + iov.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
      d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - ke * central
      cp <- central / V; cp ~ add(add.err) })
  }

  test_that("est=vae declares IOV support and theta+eta+iov stays mu-referenced", {
    expect_true(isTRUE(attr(nlmixr2est:::nlmixr2Est.vae, "iov")))
    ## the id-level eta pairs with its theta even with the extra iov term
    ui <- rxode2::assertRxUi(.vaeIovMod())
    m <- .foceiEtaThetaMap(ui)
    expect_equal(m$thetaForEta[m$etaNames == "eta.ka"], "lka")
    expect_true(is.na(m$thetaForEta[m$etaNames == "iov.ka"]))  # occasion effect, not mu-ref
  })

  test_that("est=vae trains an IOV model, holding the occasion etas fixed", {
    skip_on_cran()
    dat <- nlmixr2data::theo_md; dat$occ <- 1L; dat$occ[dat$TIME >= 144] <- 2L
    ## the IOV hook produces per-occasion fixed-variance(1) free etas
    res <- nlmixr2est:::.uiApplyIov(rxode2::assertRxUi(.vaeIovMod()), "vae", dat, vaeControl())
    .eta <- res$ui$iniDf[!is.na(res$ui$iniDf$neta1), ]
    expect_true(all(.eta$fix[grepl("^rx\\.iov", .eta$name)]))            # variance fixed
    prep <- nlmixr2est:::.vaeDataPrep(res$ui, dat)
    expect_true(all(prep$isFree[grepl("^rx\\.iov", prep$etaNames)]))     # theta forced to 0

    ctl <- vaeControl(itersBurnIn = 15L, iters = 40L, klWarmup = 12L, gammaIter = 25L,
                      nGradStep = 3L, covariateSelection = FALSE, returnVae = TRUE, seed = 1L)
    fit <- suppressMessages(suppressWarnings(nlmixr2(.vaeIovMod(), dat, est = "vae", control = ctl)))
    .iov <- grepl("^rx\\.iov", fit$prep$etaNames)
    expect_true(all(is.finite(fit$zPop)))
    expect_true(all(fit$zPop[.iov] == 0))          # occasion etas centered at 0
    expect_true(all(fit$omega[.iov] == 1))         # occasion-eta variance held at 1
    expect_true(all(fit$omega[!.iov] > 0))         # id-level omegas estimated
  })
})
