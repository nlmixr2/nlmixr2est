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

  ## ---- pinCovariates: restrict the search to model-declared covariate pairs ----

  test_that(".vaeModelCovariatePairs pairs a transformed continuous covariate", {
    ui <- rxode2::assertRxUi(.tr())               # cl.wt * log(WT/70) on cl
    dd <- as.data.frame(d); names(dd) <- toupper(names(dd))
    cov <- .vaeCovariateSearch(dd, unique(dd$ID))
    pr <- .vaeModelCovariatePairs(ui, cov$covNames, cov$covType)
    expect_equal(nrow(pr), 1L)
    expect_equal(pr$coefName, "cl.wt")
    expect_equal(pr$covName, "WT")
    expect_equal(pr$thetaName, "tcl")
    expect_equal(pr$covType, "continuous")
    expect_equal(pr$userCenter, 70)               # log(WT/70) center read from the model
    expect_true(pr$inPool)
  })

  test_that("pinCovariates=TRUE builds a 1-cell mask and zeros the training coef", {
    ui <- rxode2::assertRxUi(.tr())               # single eta (cl)
    p <- suppressWarnings(.vaeDataPrep(ui, d, vaeControl(pinCovariates = TRUE)))
    expect_true(p$pinActive)
    ## WT allowed only on the cl dim; exactly one 1 in the mask
    expect_equal(sum(p$covAllow), 1L)
    k <- match("tcl", .foceiEtaThetaMap(ui)$thetaForEta)
    expect_equal(p$covAllow[k, match("WT", p$covNames)], 1L)
    ## in-pool coefficient held at 0 for the covariate-free training decoder
    thNames <- ui$iniDf$name[!is.na(ui$iniDf$ntheta)]
    expect_equal(unname(p$th[match("cl.wt", thNames)]), 0)
    ## in-pool coefficient is NOT regressed (the search estimates it)
    expect_false("cl.wt" %in% p$regressNames)
  })

  test_that("a raw-linear effect on a continuous covariate is not pinnable (regresses)", {
    ui <- rxode2::assertRxUi(.lin())              # cl.wt * WT (linear on continuous WT)
    p <- suppressWarnings(.vaeDataPrep(ui, d, vaeControl(pinCovariates = TRUE)))
    expect_true(p$pinActive)
    expect_false(any(p$pinPairs$inPool))          # linear-on-continuous is not transferable
    expect_true("cl.wt" %in% p$regressNames)       # routed to the regress M-step
    ## mask exists and forbids every cell, so the search selects nothing
    expect_false(is.null(p$covAllow))
    expect_equal(sum(p$covAllow), 0L)
  })

  test_that("pinCovariates=FALSE turns the search off and regresses model covariates", {
    ui <- rxode2::assertRxUi(.tr())
    w <- character(0)
    p <- withCallingHandlers(
      .vaeDataPrep(ui, d, vaeControl(pinCovariates = FALSE)),
      warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
    expect_false(p$pinActive)
    expect_equal(ncol(p$covMat), 0L)               # search pool emptied
    expect_true("cl.wt" %in% p$regressNames)
    expect_true(any(grepl("pinCovariates=FALSE", w)))
  })

  test_that("a mu2 centered covariate is pinned via its nlmixrMuDerCov# column", {
    ## wt.cl*(WT/70) is a mu2 reference; the hook rewrites it to a linear
    ## nlmixrMuDerCov# column (centering carried by the mu2 data).  After the
    ## hook the covariate must be encoded LINEARLY (not log) and pinned in place.
    mu2 <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; wt.cl <- 0.1; add.err <- 0.7; eta.cl ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + wt.cl * (WT / 70) + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v; cp ~ add(add.err) })
    }
    lst <- suppressWarnings(.uiModifyForCovs(rxode2::assertRxUi(mu2()), as.data.frame(d)))
    p <- suppressWarnings(.vaeDataPrep(lst$ui, lst$data, vaeControl(pinCovariates = TRUE)))
    ## the derived column is linear (categorical), never log-encoded
    jMu <- grep("NLMIXRMUDERCOV", p$covNames)
    expect_length(jMu, 1L)
    expect_equal(p$covType[jMu], "categorical")
    ## pinned to the cl dim, in place, with no VAE-applied center (mu2 data carries it)
    expect_true(p$pinActive)
    pr <- p$pinPairs
    expect_equal(pr$coefName, "wt.cl")
    expect_true(pr$inPool)
    expect_equal(pr$userCenter, 0)
    ## searched (not regressed) and zeroed in the training theta
    expect_false("wt.cl" %in% p$regressNames)
  })

  test_that("pinCovariates=TRUE with no model covariates is a null path (full search)", {
    noCov <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.err <- 0.7; eta.cl ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v; cp ~ add(add.err) })
    }
    p <- suppressWarnings(.vaeDataPrep(rxode2::assertRxUi(noCov), d, vaeControl(pinCovariates = TRUE)))
    expect_false(p$pinActive)
    expect_null(p$covAllow)                         # no mask -> unrestricted search
    expect_equal(p$covNames, "WT")                  # WT still discovered from data
  })
})
