nmTest({
  # Older rxode2 reported `0.75*log(V1/exp(tV1))` as a mu2 covariate term with
  # covariate "-0.75"; rewriting it to `nlmixrMuDerCov*tV1` dropped
  # 0.75*log(V1), silently removing the V1 (weight, eta.V1) dependence of CL
  # for saem and the mu-referenced methods.
  test_that("a constant mu2 'covariate' is never rewritten into the model", {
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
    d <- data.frame(ID = 1, TIME = c(0, 1), AMT = c(100, 0), EVID = c(1, 0),
                    DV = c(NA, 1), WT = 80)
    lst <- .uiModifyForCovs(ui, d)
    expect_false("nlmixrMuDerCov1" %in% names(lst$data))
    .cl <- deparse1(rxode2::as.model(lst$ui)[[2]][[3]])
    expect_true(grepl("log(V1/exp(tV1))", .cl, fixed = TRUE))
    expect_false(grepl("nlmixrMuDerCov", .cl, fixed = TRUE))
  })

  test_that(".mu2CovHasVariable separates covariates from constants", {
    expect_false(.mu2CovHasVariable("-0.75"))
    expect_false(.mu2CovHasVariable("2"))
    expect_true(.mu2CovHasVariable("log(WT/70)"))
    expect_true(.mu2CovHasVariable("WT"))
  })
})
