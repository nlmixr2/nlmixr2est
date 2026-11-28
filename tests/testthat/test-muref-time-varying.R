nmTest({
  test_that("the mu-referenced regression leaves time-varying covariates out", {
    .mod <- function() {
      ini({
        lcl <- log(1)
        bwt <- 0.75
        bx <- 0.5
        eta.cl ~ 0.1
        add.sd <- 0.1
      })
      model({
        cl <- exp(lcl + bwt * WT + bx * X + eta.cl)
        cl ~ add(add.sd)
      })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.mod))
    expect_setequal(.muRefCppGroupSetup(.ui, clamp = TRUE)$muGroupCovNames, c("WT", "X"))
    # staged as time-varying (as .uiApplyMu2hook() and the imp family do)
    .nlmixrSetMuRefTimeVarying(.ui, "X")
    on.exit(.nlmixrRmMuRefTimeVarying(.ui))
    expect_equal(.muRefCovNonTv(.ui)$covariate, "WT")
    expect_equal(.muRefCppGroupSetup(.ui, clamp = TRUE)$muGroupCovNames, "WT")
  })

  test_that("mfocei estimates a time-varying mu-referenced covariate", {
    skip_on_cran()
    skip_if_not(.nlmixr2estRxode2HasTte(), "rxode2 without tte() endpoints")
    skip_if(is.null(.tteValidation()), "tte validation scenarios not installed")
    # the regression used to fit bx on the baseline X (0 for every subject)
    .fit <- .tteValidationFit("tvCov", "mfocei", list(print = 0))
    .ex <- .tteValidationExactFit("tvCov")
    expect_lt(abs(.tteValidationEst(.fit)[["bx"]] - .ex$est[["bx"]]), 0.35 * .ex$se[["bx"]])
  })
})
