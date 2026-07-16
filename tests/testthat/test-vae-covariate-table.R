## Fast regression guards for two est="vae" covariate-output bugs (no training):
##
## Bug 1 -- injected covariate-coefficient thetas (beta_<par>_<cov>) must survive
##   into the population-parameter table; they were dropped from parFixedDf when a
##   population parameter was fixed (literalFix reindex, see the slow end-to-end
##   test in test-vae-covariate.R).  Here we assert the precondition: .vaeUpdateModel
##   adds them as real thetas of the augmented ui (so they reach $theta/$cov).
##
## Bug 2 -- covariate terms must be injected FLAT (exp(theta + beta*cov + eta)),
##   not wrapped (exp((theta + beta*cov) + eta)); the wrapped form hid the exp()
##   back-transform from muRefCurEval so the mu-parameter printed on the raw log
##   scale.

nmTest({
  .theoCov <- function() {
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

  ## a minimal trained-VAE stand-in: WT selected on ka and V, not ke
  .fakeVaeFit <- function(ui) {
    prep <- .vaeDataPrep(ui, nlmixr2data::theo_sd)
    sel <- matrix(c(TRUE, FALSE, TRUE), 3, 1, dimnames = list(NULL, prep$covNames))
    beta <- matrix(c(2.5, 0, 0.5), 3, 1)
    list(prep = prep, covNames = prep$covNames, selected = sel, beta = beta,
         zPop = c(log(1.8), log(0.086), log(32)),
         omega = c(0.3, 0.03, 0.03), a = c(add.err = 0.7))
  }

  test_that("vae covariate injection is flat and preserves the exp back-transform", {
    ui <- rxode2::assertRxUi(.theoCov)
    ui2 <- suppressMessages(.vaeUpdateModel(ui, .fakeVaeFit(ui)))

    ## covariate-bearing lines are the additive mu-referenced form
    kaLine <- deparse1(ui2$lstExpr[[1]])
    vLine <- deparse1(ui2$lstExpr[[3]])
    expect_match(kaLine, "exp(lka + beta_lka_WT * log(WT/", fixed = TRUE)
    expect_match(vLine, "exp(lV + beta_lV_WT * log(WT/", fixed = TRUE)
    ## NOT the wrapped form exp((theta + beta*cov) + eta) that breaks detection
    expect_false(grepl("exp((", kaLine, fixed = TRUE))
    expect_false(grepl("exp((", vLine, fixed = TRUE))

    ## muRefCurEval still recognizes the covariate-bearing thetas as exp() so they
    ## back-transform; the non-covariate ke is unaffected
    mce <- ui2$muRefCurEval
    expect_equal(mce$curEval[mce$parameter == "lka"], "exp")
    expect_equal(mce$curEval[mce$parameter == "lV"], "exp")
    expect_equal(mce$curEval[mce$parameter == "lke"], "exp")
  })

  test_that("vae covariate coefficients become thetas of the augmented model", {
    ui <- rxode2::assertRxUi(.theoCov)
    ui2 <- suppressMessages(.vaeUpdateModel(ui, .fakeVaeFit(ui)))
    thetaNames <- ui2$iniDf$name[!is.na(ui2$iniDf$ntheta)]
    expect_true(all(c("beta_lka_WT", "beta_lV_WT") %in% thetaNames))
    ## the un-selected ke gets no coefficient
    expect_false("beta_lke_WT" %in% thetaNames)
  })

  ## A literally-fixed structural parameter (e.g. lke <- fix(...)) leaves its eta
  ## without a mu-referenced theta (thetaForEta == NA).  .vaeUpdateModel must skip
  ## that NA rather than try to set an ini() for a parameter named "NA"
  ## (previously errored: "cannot find parameter 'NA'").
  test_that("vae model update skips etas whose structural theta is fixed", {
    theoFixKe <- function() {
      ini({
        lka <- log(1.8); lke <- fix(log(0.086)); lV <- log(32)
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
    ## the model est="vae" actually sees has lke literally fixed into the RHS
    uiF <- rxode2::rxUiDecompress(
      rxode2::rxFixPop(rxode2::assertRxUi(theoFixKe), returnNull = TRUE))
    expect_true(anyNA(.foceiEtaThetaMap(uiF)$thetaForEta))   # eta.ke -> NA

    ## WT selected on ka and V, NOT on the free (fixed-theta) ke eta
    prep <- .vaeDataPrep(uiF, nlmixr2data::theo_sd)
    fit <- list(prep = prep, covNames = prep$covNames,
                selected = matrix(c(TRUE, FALSE, TRUE), 3, 1,
                                  dimnames = list(NULL, prep$covNames)),
                beta = matrix(c(2.5, 0, 0.5), 3, 1),
                zPop = prep$zPop, omega = c(0.3, 0.03, 0.03), a = c(add.err = 0.7))
    expect_error(ui2 <- suppressMessages(.vaeUpdateModel(uiF, fit)), NA)
    ## covariate coefficients still injected on the non-fixed params
    thetaNames <- ui2$iniDf$name[!is.na(ui2$iniDf$ntheta)]
    expect_true(all(c("beta_lka_WT", "beta_lV_WT") %in% thetaNames))
    ## no coefficient invented for the fixed-theta eta
    expect_false(any(grepl("beta_.*_ke|beta_NA", thetaNames)))
  })
})
