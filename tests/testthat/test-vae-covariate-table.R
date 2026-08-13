## Fast regression guards for two est="vae" covariate-output bugs (no training):
##
## Bug 1 -- injected covariate-coefficient thetas (beta.<par>.<cov>) must survive
##   into the population-parameter table; they were dropped from parFixedDf when a
##   population parameter was fixed (literalFix reindex, see the slow end-to-end
##   test in test-vae-covariate-selection.R).  Here we assert the precondition: .vaeUpdateModel
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

  ## a minimal trained-VAE stand-in: WT selected on ka and V, not ke.  `shape`
  ## picks which of WT's shape-family columns the search "selected".
  .fakeVaeFit <- function(ui, shape = "power", control = vaeControl()) {
    prep <- .vaeDataPrep(ui, nlmixr2data::theo_sd, control)
    .j <- match(shape, prep$covShape)
    sel <- matrix(FALSE, 3, length(prep$covNames),
                  dimnames = list(NULL, prep$covNames))
    beta <- matrix(0, 3, length(prep$covNames))
    sel[c(1, 3), .j] <- TRUE
    beta[c(1, 3), .j] <- c(2.5, 0.5)
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
    expect_match(kaLine, "exp(lka + beta.lka.WT.power * log(WT/", fixed = TRUE)
    expect_match(vLine, "exp(lV + beta.lV.WT.power * log(WT/", fixed = TRUE)
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
    expect_true(all(c("beta.lka.WT.power", "beta.lV.WT.power") %in% thetaNames))
    ## the un-selected ke gets no coefficient
    expect_false("beta.lke.WT.power" %in% thetaNames)
  })

  test_that("the selected shape family decides the written parameterization", {
    ui <- rxode2::assertRxUi(.theoCov)
    ui2 <- suppressMessages(.vaeUpdateModel(ui, .fakeVaeFit(ui, "lin")))
    kaLine <- deparse1(ui2$lstExpr[[1]])
    ## the linear family writes (WT - ctr), never the log form
    expect_match(kaLine, "beta.lka.WT.lin * (WT - ", fixed = TRUE)
    expect_false(grepl("log(WT", kaLine, fixed = TRUE))
    thetaNames <- ui2$iniDf$name[!is.na(ui2$iniDf$ntheta)]
    expect_true("beta.lka.WT.lin" %in% thetaNames)
    expect_false("beta.lka.WT.power" %in% thetaNames)
  })

  test_that("shapes= picks the parameterization within the selected family", {
    ui <- rxode2::assertRxUi(.theoCov)
    ## "center" spans the same model as "lin", so the SAME linear-family column
    ## is selected but written as the ratio form with a rescaled coefficient
    ctl <- vaeControl(shapes = c("center", "power"))
    fit <- .fakeVaeFit(ui, "center", ctl)
    ui2 <- suppressMessages(.vaeUpdateModel(ui, fit))
    kaLine <- deparse1(ui2$lstExpr[[1]])
    expect_match(kaLine, "beta.lka.WT.center * (WT/", fixed = TRUE)
    ## center's coefficient is the linear one times the centering value, and the
    ## structural theta absorbs the difference so the prediction is unchanged
    .j <- match("center", fit$prep$covShape)
    .ctr <- fit$prep$covPop[.j]
    .idf <- ui2$iniDf
    expect_equal(.idf$est[.idf$name == "beta.lka.WT.center"], 2.5 * .ctr)
    expect_equal(.idf$est[.idf$name == "lka"], log(1.8) - 2.5 * .ctr)
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
    fit <- .fakeVaeFit(uiF)
    fit$zPop <- fit$prep$zPop
    expect_error(ui2 <- suppressMessages(.vaeUpdateModel(uiF, fit)), NA)
    ## covariate coefficients still injected on the non-fixed params
    thetaNames <- ui2$iniDf$name[!is.na(ui2$iniDf$ntheta)]
    expect_true(all(c("beta.lka.WT.power", "beta.lV.WT.power") %in% thetaNames))
    ## no coefficient invented for the fixed-theta eta
    expect_false(any(grepl("beta\\..*\\.ke|beta\\.NA", thetaNames)))
  })
})

## Regressions for the fifth independent-review pass (write-back path).
nmTest({
  .theo <- function() {
    ini({ lka <- log(1.8); lke <- log(0.086); lV <- log(32)
      eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03; add.err <- 0.7 })
    model({ ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
      d/dt(depot) = -ka * depot
      d/dt(central) = ka * depot - ke * central
      cp <- central / V; cp ~ add(add.err) })
  }
  .fake <- function(ui, shape = "power", control = vaeControl()) {
    prep <- .vaeDataPrep(ui, nlmixr2data::theo_sd, control)
    .j <- match(shape, prep$covShape)
    sel <- matrix(FALSE, 3, length(prep$covNames))
    beta <- matrix(0, 3, length(prep$covNames))
    sel[1, .j] <- TRUE; beta[1, .j] <- 2.5
    list(prep = prep, covNames = prep$covNames, selected = sel, beta = beta,
         zPop = c(log(1.8), log(0.086), log(32)),
         omega = c(0.3, 0.03, 0.03), a = c(add.err = 0.7))
  }

  test_that("a theta named twice in its line gets the covariate injected once", {
    ## a textual gsub replaces EVERY occurrence, and operator precedence turns
    ## the second one into a spurious extra additive covariate effect
    dup <- function() {
      ini({ lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03; add.err <- 0.7 })
      model({ ka <- exp(lka + eta.ka) + 0 * lka
        ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot
        d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ui <- rxode2::assertRxUi(dup)
    ui2 <- suppressMessages(.vaeUpdateModel(ui, .fake(ui)))
    .line <- deparse1(ui2$lstExpr[[1]])
    expect_equal(lengths(regmatches(.line, gregexpr("beta.lka.WT.power", .line))), 1L)
    ## still the flat mu-referenced form, so the exp() back-transform survives
    expect_match(.line, "exp(lka + beta.lka.WT.power * log(WT/", fixed = TRUE)
    expect_equal(ui2$muRefCurEval$curEval[ui2$muRefCurEval$parameter == "lka"], "exp")
  })

  test_that("an intercept correction never writes a theta out of bounds", {
    ## "identity" moves beta*center into the intercept; with a bounded theta the
    ## corrected value can leave the bounds, and clamping it would change the
    ## prediction -- so the centered parameterization is written instead
    bnd <- function() {
      ini({ lka <- c(0.5, log(1.8), 0.7); lke <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03; add.err <- 0.7 })
      model({ ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot
        d/dt(central) = ka * depot - ke * central
        cp <- central / V; cp ~ add(add.err) })
    }
    ui <- rxode2::assertRxUi(bnd)
    ctl <- vaeControl(shapes = c("identity", "power"))
    ui2 <- suppressMessages(.vaeUpdateModel(ui, .fake(ui, "identity", ctl)))
    .idf <- ui2$iniDf
    .est <- .idf$est[.idf$name == "lka"]
    ## the written estimate respects the declared bounds
    expect_gt(.est, 0.5)
    expect_lt(.est, 0.7)
    ## and it stayed exact by writing the centered form rather than clamping
    expect_true(any(grepl("beta.lka.WT.lin", .idf$name)))
    expect_equal(.est, log(1.8))
  })

  test_that("generated coefficient names stay distinct", {
    expect_equal(.vaeUniqueName("beta.lka.WT.log", character(0)), "beta.lka.WT.log")
    expect_equal(.vaeUniqueName("beta.lka.WT.log", "beta.lka.WT.log"),
                 "beta.lka.WT.log2")
    expect_equal(.vaeUniqueName("beta.lka.WT.log",
                                c("beta.lka.WT.log", "beta.lka.WT.log2")),
                 "beta.lka.WT.log3")
  })
})
