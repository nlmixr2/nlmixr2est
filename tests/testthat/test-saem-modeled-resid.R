nmTest({
  # rxUse references every ini() parameter so each variant parses
  .modeledResidUi <- function(errLine, body = "", ini = "") {
    .txt <- sprintf("function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        add.sd <- 0.7; prop.sd <- 0.1; pw <- 1; lam <- 0.5; cov.sd <- 0.001
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.sd ~ 0.1
        %s
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        cp <- linCmt()
        rxUse1 <- add.sd; rxUse2 <- prop.sd; rxUse3 <- pw; rxUse4 <- lam
        rxUse5 <- cov.sd; rxUse6 <- exp(eta.sd)
        %s
        %s
      })
    }", ini, body, errLine)
    rxode2::rxode2(eval(parse(text = .txt)))
  }
  .hook <- function(ui) {
    suppressWarnings(.preProcessSaemModeledResid(ui, "saem", NULL, NULL)$ui)
  }
  .lines <- function(ui) vapply(ui$lstExpr, deparse1, character(1))

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
    expect_null(.preProcessSaemModeledResid(.ui, "saem", NULL, NULL))
  })

  test_that("the saem hook rewrites a modeled residual to its + dnorm() twin", {
    nlmixr2global$nlmixr2EstEnv$nlmixrPureInputUi <- NULL
    on.exit(nlmixr2global$nlmixr2EstEnv$nlmixrPureInputUi <- NULL, add = TRUE)
    .ui <- .modeledResidUi("cp ~ add(a)", "a <- add.sd * exp(eta.sd)")
    expect_null(.preProcessSaemModeledResid(.ui, "focei", NULL, NULL))
    expect_warning(.new <- .preProcessSaemModeledResid(.ui, "saem", NULL, NULL)$ui,
                   "modeled residual error for 'cp'")
    expect_equal(as.character(.new$predDf$distribution), "dnorm")
    expect_equal(.new$saemResMod, 0L)
    expect_true("a <- add.sd * exp(eta.sd)" %in% .lines(.new))
    expect_match(paste(deparse(.new$saemModel0), collapse = "\n"), "llikNorm", fixed = TRUE)
    # the model as written is kept for the reported fit
    expect_equal(nlmixr2global$nlmixr2EstEnv$nlmixrPureInputUi$lstExpr, .ui$lstExpr)
    nlmixr2global$nlmixr2EstEnv$nlmixrPureInputUi <- NULL
    .twin <- .hook(.modeledResidUi("cp ~ add(a) + dnorm()", "a <- add.sd * exp(eta.sd)"))
    expect_equal(.new$saemModel0, .twin$saemModel0)
  })

  test_that("a covariate on a residual theta without an eta stays in the model", {
    .ui <- .modeledResidUi("cp ~ add(a)", "a <- add.sd + WT * cov.sd")
    expect_equal(.ui$muRefCovariateDataFrame$covariateParameter, "cov.sd")
    .new <- .hook(.ui)
    expect_equal(nrow(.new$saemMuRefCovariateDataFrame), 0L)
    expect_true("WT" %in% .new$saemInPars)
    expect_true(all(c("rxBoundedTr.add.sd", "rxBoundedTr.cov.sd") %in% .new$saemParamsToEstimate))
    expect_true("a <- add.sd + WT * cov.sd" %in% .lines(.new))
  })

  test_that("a modeled residual with a lambda transform is refused", {
    .ui <- .modeledResidUi("cp ~ add(add.sd) + boxCox(l2)", "l2 <- lam + WT * cov.sd")
    expect_error(.preProcessSaemModeledResid(.ui, "saem", NULL, NULL),
                 "boxCox()/yeoJohnson()", fixed = TRUE)
    .ui <- .modeledResidUi("cp ~ add(a) + yeoJohnson(lam)", "a <- add.sd * exp(eta.sd)")
    expect_error(.preProcessSaemModeledResid(.ui, "saem", NULL, NULL),
                 "boxCox()/yeoJohnson()", fixed = TRUE)
  })

  test_that("propF()/powF() prediction variables do not promote the endpoint", {
    .ui <- .modeledResidUi("cp ~ propF(prop.sd, f2)", "f2 <- 1")
    expect_equal(.saemModeledResidualCond(.ui), character(0))
    .ui <- .modeledResidUi("cp ~ powF(prop.sd, pw, f2)", "f2 <- 1")
    expect_equal(.saemModeledResidualCond(.ui), character(0))
  })

  test_that("the temporary eta scale follows the likelihood argument's range", {
    # a theta that IS the argument takes the argument's range
    .s <- .saemPseudoEtaThetas(.modeledResidUi("cp ~ add(add.sd) + dnorm()"))
    expect_equal(.s$theta, "add.sd")
    expect_equal(c(.s$lower, .s$upper), c(0, Inf))
    .s <- .saemPseudoEtaThetas(.hook(.modeledResidUi("cp ~ add(a)", "a <- add.sd")))
    expect_equal(nrow(.s), 0L) # already carries its temporary eta
    .new <- .hook(.modeledResidUi("cp ~ add(a)", "a <- add.sd"))
    expect_true("add.sd <- exp(rxBoundedTr.add.sd + rx.eta.add.sd)" %in% .lines(.new))
    expect_equal(.new$iniDf$est[.new$iniDf$name == "rxBoundedTr.add.sd"], log(0.7))
    expect_true(any(.new$muRefDataFrame$theta == "rxBoundedTr.add.sd" &
                      .new$muRefDataFrame$eta == "rx.eta.add.sd"))
    expect_equal(.new$boundedTransforms[[1]]$type, "lower_exp")
    # a theta inside an expression uses its own bounds: unbounded -> additive
    .new <- .hook(.modeledResidUi("cp ~ add(a)", "a <- add.sd * exp(eta.sd)"))
    expect_true("add.sd <- rxBoundedTr.add.sd + rx.eta.add.sd" %in% .lines(.new))
    expect_equal(.new$boundedTransforms[[1]]$type, "identity")
    # ... and a bounded theta inside an expression uses those bounds
    .new <- .hook(.modeledResidUi("cp ~ add(a)", "a <- add.sd * exp(eta.sd) + WT * cov2",
                                  ini = "cov2 <- c(0, 0.01, 1)"))
    expect_true("cov2 <- expit(rxBoundedTr.cov2 + rx.eta.cov2, 0, 1)" %in% .lines(.new))
    # t() degrees of freedom and a binomial probability
    .ui <- .modeledResidUi("cp ~ add(add.sd) + dt(nu)", ini = "nu <- 5")
    .s <- .saemPseudoEtaThetas(.ui)
    expect_equal(c(.s$lower[.s$theta == "nu"], .s$upper[.s$theta == "nu"]), c(0, Inf))
    .ub <- rxode2::rxode2(function() {
      ini({ tp <- 0; p <- 0.3; eta.p ~ 0.1 })
      model({ q <- expit(tp + eta.p); dv ~ dbinom(10, p) })
    })
    .s <- .saemPseudoEtaThetas(.ub)
    expect_equal(.s$theta, "p")
    expect_equal(c(.s$lower, .s$upper), c(0, 1))
  })

  test_that("temporary eta lines stay mu-referenced for every range", {
    .ui <- .modeledResidUi("cp ~ add(add.sd) + dnorm()")
    for (.r in list(c(0, Inf), c(0.1, Inf), c(-Inf, 5), c(0.1, 3), c(-Inf, Inf))) {
      .new <- .saemAddPseudoEtas(.ui, data.frame(theta = "add.sd", lower = .r[1], upper = .r[2]))
      expect_true(any(.new$muRefDataFrame$eta == "rx.eta.add.sd"), info = paste(.r, collapse = ","))
      expect_false("rx.eta.add.sd" %in% .new$nonMuEtas, info = paste(.r, collapse = ","))
    }
    .new <- .saemAddPseudoEtas(.ui, data.frame(theta = "add.sd", lower = 0.1, upper = Inf))
    expect_true(all(c("rx.l.add.sd <- exp(rxBoundedTr.add.sd + rx.eta.add.sd)",
                      "add.sd <- 0.1 + rx.l.add.sd") %in% .lines(.new)))
  })

  test_that("temporary etas skip mu-referenced, fixed, count and prediction thetas", {
    expect_equal(nrow(.saemPseudoEtaThetas(.hook(
      .modeledResidUi("cp ~ add(a)", "a <- exp(lam + eta.sd)")))), 0L)
    .s <- .saemPseudoEtaThetas(.modeledResidUi("cp ~ add(a) + dnorm()", "a <- add.sd * cp"))
    expect_equal(.s$theta, "add.sd")
    .uf <- rxode2::rxode2(function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- fixed(0.7); eta.ka ~ 0.6; eta.cl ~ 0.3 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
              cp <- linCmt(); cp ~ add(add.sd) + dnorm() })
    })
    expect_equal(nrow(.saemPseudoEtaThetas(.uf)), 0L)
    .ub <- rxode2::rxode2(function() {
      ini({ tp <- 0; n <- 10; p <- 0.3; eta.p ~ 0.1 })
      model({ q <- expit(tp + eta.p); dv ~ dbinom(n, p) })
    })
    expect_equal(.saemPseudoEtaThetas(.ub)$theta, "p")
    expect_equal(nrow(.saemPseudoEtaThetas(.modeledResidUi("cp ~ add(add.sd)"))), 0L)
  })

  test_that("ll() thetas informing the likelihood get temporary etas", {
    .ul <- rxode2::rxode2(function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45; tka.wt <- 0.01; lsd <- log(0.7)
        eta.ka ~ 0.6; eta.cl ~ 0.3
      })
      model({
        ka <- exp(tka + tka.wt * WT + eta.ka)
        cl <- exp(tcl + eta.cl); v <- exp(tv)
        sd <- exp(lsd)
        cp <- linCmt()
        ll(err) ~ -log(sd) - 0.5 * log(2 * pi) - 0.5 * ((DV - cp) / sd)^2
      })
    })
    .s <- .saemPseudoEtaThetas(.ul)
    # the scale and the eta-less structural theta; not tka (has an eta) or its covariate
    expect_setequal(.s$theta, c("lsd", "tv"))
    expect_true(all(is.infinite(c(.s$lower, .s$upper))))
    .new <- .hook(.ul)
    expect_true(all(c("rx.eta.lsd", "rx.eta.tv") %in% .new$iniDf$name))
    expect_true("tka.wt" %in% .new$iniDf$name)
    expect_length(.new$nonMuEtas, 0L)
  })

  test_that("each modeled endpoint's residual theta gets a temporary eta", {
    .u <- rxode2::rxode2(function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add1 <- 0.1; add2 <- 0.2
            eta.ka ~ 0.6; eta.cl ~ 0.3; eta1 ~ 0.1; eta2 ~ 0.1 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
              cp <- linCmt(); cp2 <- 2 * cp
              a1 <- add1 * exp(eta1)
              cp ~ add(a1)
              a2 <- add2 * exp(eta2)
              cp2 ~ add(a2) })
    })
    .new <- .hook(.u)
    expect_equal(as.character(.new$predDf$distribution), c("dnorm", "dnorm"))
    expect_equal(sum(grepl("~.*dnorm\\(\\)", .lines(.new))), 2L)
    # the first error line must not read as an assignment that hides add1
    expect_true(all(c("rx.eta.add1", "rx.eta.add2") %in% .new$iniDf$name))
  })

  test_that("temporary-eta transforms survive the bounded-transform hook", {
    .u <- rxode2::rxode2(function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; t1 <- c(0, 1, 2); add.sd <- 0.1
            eta.ka ~ 0.6; eta.cl ~ 0.3; eta.sd ~ 0.1 })
      model({ ka <- exp(tka + eta.ka) * t1; cl <- exp(tcl + eta.cl); v <- exp(tv)
              cp <- linCmt(); a <- add.sd * exp(eta.sd); cp ~ add(a) })
    })
    .new <- .hook(.u)
    .stash <- .new$boundedTransforms
    .bt <- suppressWarnings(.preProcessBoundedTransform(.new, "saem", NULL, saemControl())$ui)
    .names <- function(ui) vapply(ui$boundedTransforms, function(tr) tr$name, character(1))
    # the bounded-transform hook keeps only the user's bounded theta ...
    expect_false("add.sd" %in% .names(.bt))
    # ... so saem adds its own back without losing the user's
    .bt <- .saemRestorePseudoTransforms(.bt, .stash)
    expect_setequal(.names(.bt), c("t1", "add.sd"))
    .bt <- .saemRestorePseudoTransforms(.bt, .stash)
    expect_length(.bt$boundedTransforms, 2L)
    # a ui that already carries the component (set directly, as nlmixr2Est0 does) and is
    # then compressed must still accept the merge
    .d <- rxode2::rxUiDecompress(suppressWarnings(
      .preProcessBoundedTransform(.new, "saem", NULL, saemControl())$ui))
    .d$boundedTransforms <- .d$boundedTransforms
    .d <- .saemRestorePseudoTransforms(rxode2::rxUiCompress(.d), .stash)
    expect_setequal(.names(.d), c("t1", "add.sd"))
  })

  test_that("a bounded structural theta and a temporary eta are both back-transformed", {
    .m <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; t1 <- c(0, 1, 2); add.sd <- 0.7
            eta.ka ~ 0.6; eta.cl ~ 0.3; eta.sd ~ 0.1 })
      model({ ka <- exp(tka + eta.ka) * t1; cl <- exp(tcl + eta.cl); v <- exp(tv)
              cp <- linCmt(); a <- add.sd * exp(eta.sd); cp ~ add(a) })
    }
    .f <- suppressWarnings(.nlmixr(.m, nlmixr2data::theo_sd, est = "saem",
      control = saemControl(nBurn = 10, nEm = 10, seed = 42L, print = 0L, covMethod = "",
                            calcTables = FALSE)))
    expect_true(all(c("t1", "add.sd") %in% names(fixef(.f))))
    expect_false(any(grepl("^rx", names(fixef(.f)))))
    expect_true(fixef(.f)[["t1"]] > 0 && fixef(.f)[["t1"]] < 2)
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
    fD <- .nlmixr(mDnorm, nlmixr2data::theo_sd, est = "saem", control = ctl)
    expect_true(any(grepl("modeled residual error for 'cp'", fM$runInfo, fixed = TRUE)))
    expect_true(any(grepl("temporary eta for eta-less likelihood theta(s): add.sd",
                          fM$runInfo, fixed = TRUE)))
    expect_equal(fixef(fM), fixef(fD))
    expect_equal(fM$omega, fD$omega)
    expect_equal(fM$objf, fD$objf)
    for (.f in list(fM, fD)) {
      # the temporary eta is gone and the model is reported as written
      expect_false(any(grepl("^rx[.](eta|l)[.]|^rxBoundedTr", c(names(.f), names(.f$eta),
                                                               rownames(.f$omega),
                                                               names(fixef(.f))))))
      expect_true("add.sd" %in% names(fixef(.f)))
      expect_false(any(grepl("^rx", .f$finalUi$iniDf$name)))
      # the table reports the prediction, not the log-density (#1084)
      .d <- as.data.frame(.f)
      expect_equal(.d$IPRED, .d$cp)
      expect_equal(.d$IRES, .d$DV - .d$IPRED)
      expect_equal(.d$IWRES, .d$IRES / .d$a)
    }
    expect_equal(fM$finalUi$lstExpr, rxode2::rxode2(mMod)$lstExpr)
  })
})
