nmTest({
  # rxUse references every ini() parameter so each variant parses
  .modeledResidUi <- function(errLine, body = "") {
    .txt <- sprintf("function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        add.sd <- 0.7; prop.sd <- 0.1; pw <- 1; lam <- 0.5; cov.sd <- 0.001
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.sd ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        cp <- linCmt()
        rxUse1 <- add.sd; rxUse2 <- prop.sd; rxUse3 <- pw; rxUse4 <- lam
        rxUse5 <- cov.sd; rxUse6 <- exp(eta.sd)
        %s
        %s
      })
    }", body, errLine)
    rxode2::rxode2(eval(parse(text = .txt)))
  }

  test_that("modeled residual components are detected per endpoint", {
    .ui <- .modeledResidUi("cp ~ add(a)", "a <- add.sd * exp(eta.sd)")
    expect_equal(.saemModeledResidualCond(.ui), "cp")
    .ui <- .modeledResidUi("cp ~ add(a)", "a <- add.sd + WT * cov.sd")
    expect_equal(.saemModeledResidualCond(.ui), "cp")
    .ui <- .modeledResidUi("cp ~ add(add.sd) + prop(b)", "b <- prop.sd * exp(eta.sd)")
    expect_equal(.saemModeledResidualCond(.ui), "cp")
    .ui <- .modeledResidUi("cp ~ add(add.sd) + pow(prop.sd, c2)", "c2 <- pw + cov.sd * WT")
    expect_equal(.saemModeledResidualCond(.ui), "cp")
    .ui <- .modeledResidUi("cp ~ add(add.sd) + boxCox(l2)", "l2 <- lam + WT * cov.sd")
    expect_equal(.saemModeledResidualCond(.ui), "cp")
  })

  test_that("plain theta residual components are not promoted", {
    .ui <- .modeledResidUi("cp ~ add(add.sd) + prop(prop.sd)")
    expect_equal(.saemModeledResidualCond(.ui), character(0))
    expect_false(.saemGeneralLik(.ui))
    expect_equal(.ui$saemResMod, 4L)
  })

  test_that("the saem hook rewrites a modeled residual to its + dnorm() twin", {
    .ui <- .modeledResidUi("cp ~ add(a)", "a <- add.sd * exp(eta.sd)")
    expect_null(.preProcessSaemModeledResid(.ui, "focei", NULL, NULL))
    expect_null(.preProcessSaemModeledResid(
      .modeledResidUi("cp ~ add(add.sd) + prop(prop.sd)"), "saem", NULL, NULL))
    expect_warning(.new <- .preProcessSaemModeledResid(.ui, "saem", NULL, NULL)$ui,
                   "modeled residual error for 'cp'")
    expect_equal(as.character(.new$predDf$distribution), "dnorm")
    expect_equal(.new$saemResMod, 0L)
    .txt <- paste(deparse(.new$saemModel0), collapse = "\n")
    expect_match(.txt, "llikNorm", fixed = TRUE)
    expect_match(.txt, "a <- add.sd * exp(eta.sd)", fixed = TRUE)
    .uiD <- .modeledResidUi("cp ~ add(a) + dnorm()", "a <- add.sd * exp(eta.sd)")
    expect_equal(.new$saemModel0, .uiD$saemModel0)
  })

  test_that("a modeled residual with a lambda transform is refused", {
    .ui <- .modeledResidUi("cp ~ add(add.sd) + boxCox(l2)", "l2 <- lam + WT * cov.sd")
    expect_error(.preProcessSaemModeledResid(.ui, "saem", NULL, NULL),
                 "boxCox()/yeoJohnson()", fixed = TRUE)
    .ui <- .modeledResidUi("cp ~ add(a) + yeoJohnson(lam)", "a <- add.sd * exp(eta.sd)")
    expect_error(.preProcessSaemModeledResid(.ui, "saem", NULL, NULL),
                 "boxCox()/yeoJohnson()", fixed = TRUE)
  })

  test_that("dnorm() is inserted before a | condition", {
    expect_equal(.saemAddDnormToErrLine(quote(cp ~ add(a) | cp)),
                 quote(cp ~ add(a) + dnorm() | cp))
    expect_equal(.saemAddDnormToErrLine(quote(cp ~ add(add.sd) + prop(b))),
                 quote(cp ~ add(add.sd) + prop(b) + dnorm()))
  })

  mMod <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.sd ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      cp <- linCmt()
      a <- add.sd * exp(eta.sd)
      cp ~ add(a)
    })
  }
  mDnorm <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.sd ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      cp <- linCmt()
      a <- add.sd * exp(eta.sd)
      cp ~ add(a) + dnorm()
    })
  }

  test_that("saem fits a modeled residual exactly as its + dnorm() twin", {
    ctl <- saemControl(nBurn = 20, nEm = 20, seed = 42L, print = 0L, covMethod = "")
    fM <- .nlmixr(mMod, nlmixr2data::theo_sd, est = "saem", control = ctl)
    fD <- suppressWarnings(.nlmixr(mDnorm, nlmixr2data::theo_sd, est = "saem", control = ctl))
    expect_equal(fM$ui$saemResMod, 0L)
    expect_true(any(grepl("modeled residual error for 'cp'", fM$runInfo, fixed = TRUE)))
    expect_equal(fixef(fM), fixef(fD))
    expect_equal(fM$omega, fD$omega)
    expect_equal(fM$objf, fD$objf)
    # the table reports the prediction, not the log-density (#1084)
    for (.f in list(fM, fD)) {
      .d <- as.data.frame(.f)
      expect_equal(.d$IPRED, .d$cp)
      expect_equal(.d$IRES, .d$DV - .d$IPRED)
      expect_equal(.d$IWRES, .d$IRES / .d$a)
    }
  })
})
