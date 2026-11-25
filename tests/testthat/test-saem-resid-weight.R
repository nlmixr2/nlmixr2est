nmTest({
  .weightUi <- function(errLine, body = "") {
    .txt <- sprintf(
      "function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        add.sd <- 0.7; cov.sd <- 0.1
        eta.ka ~ 0.6; eta.cl ~ 0.3
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        cp <- linCmt()
        rxUse1 <- cov.sd
        %s
        %s
      })
    }",
      body,
      errLine
    )
    rxode2::rxode2(eval(parse(text = .txt)))
  }
  .specs <- function(ui) {
    .saemResidWeightSpecs(ui, .saemModeledResidualCond(ui))
  }

  test_that("a theta times a known weight is factored out (#1463)", {
    expect_equal(.saemFactorTheta(quote(add.sd * exp(-time / 8)), "add.sd"), quote(exp(-time / 8)))
    expect_equal(.saemFactorTheta(quote(add.sd * WT / 70), "add.sd"), quote(WT / 70))
    expect_equal(.saemFactorTheta(quote(exp(-time) * (add.sd)), "add.sd"), quote(exp(-time)))
    expect_null(.saemFactorTheta(quote(WT / add.sd), "add.sd"))
    expect_null(.saemFactorTheta(quote(add.sd + time), "add.sd"))
    expect_null(.saemFactorTheta(quote(add.sd * exp(-time / add.sd)), "add.sd"))
  })

  test_that("only add(theta * w(time, covariates)) qualifies (#1463)", {
    .s <- .specs(.weightUi("cp ~ add(add.sd * exp(-time / 8))"))
    expect_length(.s, 1L)
    expect_equal(.s[[1]]$theta, "add.sd")
    expect_equal(.s[[1]]$cond, "cp")
    expect_length(.specs(.weightUi("cp ~ add(a)", "a <- add.sd * WT / 70")), 1L)
    # an eta, a second theta or an unknown function stays on the dnorm() path
    expect_length(.specs(.weightUi("cp ~ add(a)", "a <- add.sd * exp(eta.cl)")), 0L)
    expect_length(.specs(.weightUi("cp ~ add(a)", "a <- add.sd * exp(-time * cov.sd)")), 0L)
    expect_length(.specs(.weightUi("cp ~ add(a)", "a <- add.sd * expit(time)")), 0L)
    expect_length(.specs(.weightUi("cp ~ add(a)", "a <- add.sd + time")), 0L)
    # the weight may not depend on a model variable
    expect_length(.specs(.weightUi("cp ~ add(a)", "a <- add.sd * cl")), 0L)
    expect_length(.specs(.weightUi("cp ~ add(add.sd) + prop(cov.sd)")), 0L)
  })

  test_that("the hook fits a weighted residual as a scaled additive one (#1463)", {
    nlmixr2global$nlmixr2EstEnv$nlmixrPureInputUi <- NULL
    on.exit(nlmixr2global$nlmixr2EstEnv$nlmixrPureInputUi <- NULL, add = TRUE)
    .ui <- .weightUi("cp ~ add(add.sd * exp(-time / 8))")
    expect_warning(.ret <- .preProcessSaemModeledResid(.ui, "saem", NULL, NULL), NA)
    .new <- .ret$ui
    expect_false(.saemGeneralLik(.new))
    expect_equal(.new$saemResNames, "add.sd")
    expect_equal(.new$predDf$cond, "cp")
    expect_equal(paste(.new$predDf$var), "rx.cp.wpred")
    expect_false(any(grepl("rx.eta.", .new$iniDf$name, fixed = TRUE)))
    .specs <- nlmixr2global$nlmixr2EstEnv$saemResidWeights
    expect_length(.specs, 1L)
    .back <- .saemResidWeightRestoreUi(.new, .specs)
    expect_equal(
      vapply(.back$lstExpr, deparse1, character(1)),
      vapply(.ui$lstExpr, deparse1, character(1))
    )
    expect_true(is.na(.back$iniDf$err[.back$iniDf$name == "add.sd"]))
  })

  test_that("the weighted endpoint's DV and LIMIT are scaled (#1463)", {
    .d <- data.frame(
      ID = 1,
      TIME = c(0, 1, 2, 3),
      EVID = c(1, 0, 0, 0),
      CMT = c(1, 2, 2, 2),
      DV = c(NA, 4, 4, 4),
      LIMIT = c(NA, -Inf, 2, 2),
      WT = 140
    )
    .s <- list(list(cond = "cp", cmt = 2, w = quote(WT / 70 * exp(-time))))
    .r <- .saemResidWeightData(.d, .s, 1L)
    .w <- 2 * exp(-(1:3))
    expect_equal(.r$data$DV, c(NA, 4 / .w))
    expect_equal(.r$data$LIMIT, c(NA, -Inf, 2 / .w[2:3]))
    expect_equal(.r$logW, sum(log(.w)))
    expect_identical(.saemResidWeightData(.d, list(), 1L)$data, .d)
    .s[[1]]$w <- quote(time - 2)
    expect_error(.saemResidWeightData(.d, .s, 1L), "must be positive")
  })
})
