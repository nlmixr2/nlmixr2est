## With covariateSelection=FALSE the covariate coefficients WRITTEN in the model
## (linear beta*WT and transformed beta*log(WT/70)) must be estimated by the
## regress M-step regardless of nonMuTheta -- previously they were held fixed
## under "none" and errored under "fix"/"eta".  Covariate identity is read from
## the shared muRefCovariateDataFrame/allCovs (never mutated).  Fast: no fit; the
## end-to-end "coefficient moves" checks live in the slow test-vae-covariate.R.

nmTest({
  ## linear (rxode2 records it in muRefCovariateDataFrame) covariate on cl
  .lin <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; cl.wt <- 0.01; add.err <- 0.7; eta.cl ~ 0.1 })
    model({ ka <- exp(tka); cl <- exp(tcl + cl.wt * WT + eta.cl); v <- exp(tv)
      d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v; cp ~ add(add.err) })
  }
  ## transformed (NOT recognized as a mu-ref covariate) effect on cl
  .tr <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; cl.wt <- 0.1; add.err <- 0.7; eta.cl ~ 0.1 })
    model({ ka <- exp(tka); cl <- exp(tcl + cl.wt * log(WT / 70) + eta.cl); v <- exp(tv)
      d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v; cp ~ add(add.err) })
  }
  d <- nlmixr2data::theo_sd

  test_that(".vaeCovariateCoefThetas detects linear and transformed coefficients", {
    expect_equal(.vaeCovariateCoefThetas(rxode2::assertRxUi(.lin())), "cl.wt")
    expect_equal(.vaeCovariateCoefThetas(rxode2::assertRxUi(.tr())), "cl.wt")
  })

  test_that(".vaeCovariateCoefThetas drops a user-fixed coefficient", {
    fx <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; cl.wt <- fix(0.01); add.err <- 0.7; eta.cl ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + cl.wt * WT + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v; cp ~ add(add.err) })
    }
    expect_equal(.vaeCovariateCoefThetas(rxode2::assertRxUi(fx())), character(0))
  })

  test_that(".vaeNonMuThetas excludes covariate coefficients (they are not structural)", {
    ## cl.wt must NOT be offered to nonMuTheta eta/fix injection
    expect_false("cl.wt" %in% .vaeNonMuThetas(rxode2::assertRxUi(.lin())))
    expect_false("cl.wt" %in% .vaeNonMuThetas(rxode2::assertRxUi(.tr())))
  })

  test_that("covariateSelection=FALSE regresses the coefficient in every nonMuTheta mode", {
    for (m in c("regress", "none", "fix", "eta")) {
      p <- suppressWarnings(.vaeDataPrep(rxode2::assertRxUi(.tr()), d,
                                         vaeControl(covariateSelection = FALSE, nonMuTheta = m)))
      expect_true("cl.wt" %in% p$regressNames, info = m)
      i <- match("cl.wt", p$regressNames)
      expect_true(is.finite(p$regressLower[i]) && is.finite(p$regressUpper[i]), info = m)
    }
  })

  test_that("covariateSelection=TRUE does NOT force the coefficient into the regress set", {
    ## default nonMuTheta='regress' still picks up genuine non-mu structural thetas
    ## (tka, tv) but the coefficient is left to the selection machinery
    p <- suppressWarnings(.vaeDataPrep(rxode2::assertRxUi(.tr()), d,
                                       vaeControl(covariateSelection = TRUE, nonMuTheta = "none")))
    expect_false("cl.wt" %in% p$regressNames)
  })

  test_that("a fixed coefficient is not regressed", {
    fx <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; cl.wt <- fix(0.01); add.err <- 0.7; eta.cl ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + cl.wt * WT + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v; cp ~ add(add.err) })
    }
    p <- suppressWarnings(.vaeDataPrep(rxode2::assertRxUi(fx()), d,
                                       vaeControl(covariateSelection = FALSE, nonMuTheta = "none")))
    expect_false("cl.wt" %in% p$regressNames)
  })

  test_that("the unbounded fallback bound is scale-aware (tighter for a raw covariate)", {
    ## raw WT (~O(70)) => bound ~ 1/max|WT|, much tighter than the log-scale default
    pLin <- suppressWarnings(.vaeDataPrep(rxode2::assertRxUi(.lin()), d,
                                          vaeControl(covariateSelection = FALSE, nonMuTheta = "none")))
    pTr <- suppressWarnings(.vaeDataPrep(rxode2::assertRxUi(.tr()), d,
                                         vaeControl(covariateSelection = FALSE, nonMuTheta = "none")))
    bLin <- pLin$regressUpper[match("cl.wt", pLin$regressNames)]
    bTr <- pTr$regressUpper[match("cl.wt", pTr$regressNames)]
    expect_lt(bLin, bTr)
    expect_equal(bLin, .vaeCovCoefEffect / max(abs(d$WT)))
  })

  test_that("the fit surfaces a 'estimating covariate coef(s)' note and never errors", {
    for (m in c("none", "fix", "eta", "regress")) {
      w <- character(0)
      withCallingHandlers(
        .preProcessVaeNonMuTheta(rxode2::assertRxUi(.tr()), "vae", d,
                                 vaeControl(covariateSelection = FALSE, nonMuTheta = m)),
        warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
      expect_true(any(grepl("estimating covariate coef", w)), info = m)
    }
  })
})
