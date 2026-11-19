nmTest({
  # rxode2 releases without nlmixr2/rxode2#1756 report `0.75*log(V1/exp(tV1))`
  # as a mu2 covariate term with covariate "-0.75".  Rewriting it to
  # `nlmixrMuDerCov*tV1` dropped 0.75*log(V1), silently removing the V1
  # (weight, eta.V1) dependence of CL for saem and the mu-referenced methods.
  .constantMu2Ui <- function() {
    f <- function() {
      ini({
        tV1 <- log(40); tCL <- log(5)
        eta.CL ~ 0.1; eta.V1 ~ 0.1
        add.sd <- 1
      })
      model({
        V1 <- exp(tV1 + log(WT/70) + eta.V1)
        CL <- exp(tCL + 0.75*log(V1/exp(tV1)) + eta.CL)
        cp <- linCmt(CL, V1)
        cp ~ add(add.sd)
      })
    }
    ui <- rxode2::rxUiDecompress(rxode2::assertRxUi(f))
    # inject what the older rxode2 produced, independent of the installed version
    assign("mu2RefCovariateReplaceDataFrame",
           data.frame(theta = "tCL", covariate = "-0.75", covariateParameter = "tV1",
                      modelExpression = "0.75 * log(V1/exp(tV1))"),
           envir = ui)
    ui
  }

  test_that("a constant mu2 'covariate' is never rewritten into the model", {
    ui <- .constantMu2Ui()
    d <- data.frame(ID = 1, TIME = c(0, 1), AMT = c(100, 0), EVID = c(1, 0),
                    DV = c(NA, 1), WT = 80)
    expect_warning(lst <- .uiModifyForCovs(ui, d), "constant covariate")
    expect_false("nlmixrMuDerCov1" %in% names(lst$data))
    .cl <- deparse1(rxode2::as.model(lst$ui)[[2]][[3]])
    expect_true(grepl("log(V1/exp(tV1))", .cl, fixed = TRUE))
    expect_false(grepl("nlmixrMuDerCov", .cl, fixed = TRUE))
  })

  test_that("the pre-processing hook drops constant mu2 rows for every reader", {
    ui <- .constantMu2Ui()
    expect_warning(ret <- .nlmixr0preProcessMu2ConstantCov(ui, "saem", NULL, NULL),
                   "0.75 \\* log\\(V1/exp\\(tV1\\)\\)")
    expect_named(ret, "ui")
    expect_equal(nrow(ret$ui$mu2RefCovariateReplaceDataFrame), 0L)
    # the user's ui is not modified
    expect_equal(nrow(ui$mu2RefCovariateReplaceDataFrame), 1L)
    # it runs before the mu2 rewrite
    .hooks <- preProcessHooks()
    expect_lt(match(".nlmixr0preProcessMu2ConstantCov", .hooks), match(".uiApplyMu2", .hooks))
  })

  test_that("the hook leaves genuine mu2 covariates alone", {
    f <- function() {
      ini({ tCL <- log(5); tV <- log(40); cl.wt <- 0.75; eta.CL ~ 0.1; add.sd <- 1 })
      model({
        CL <- exp(tCL + cl.wt*log(WT/70) + eta.CL)
        V <- exp(tV)
        cp <- linCmt(CL, V)
        cp ~ add(add.sd)
      })
    }
    ui <- rxode2::assertRxUi(f)
    expect_null(.nlmixr0preProcessMu2ConstantCov(ui, "saem", NULL, NULL))
  })

  test_that(".mu2CovHasVariable separates covariates from constants", {
    expect_false(.mu2CovHasVariable("-0.75"))
    expect_false(.mu2CovHasVariable("2"))
    expect_true(.mu2CovHasVariable("log(WT/70)"))
    expect_true(.mu2CovHasVariable("WT"))
  })
})
