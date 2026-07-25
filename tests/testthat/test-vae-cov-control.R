nmTest({
  test_that("vaeControl exposes the covariate shape settings", {
    .c <- vaeControl()
    expect_equal(.c$shapes,
                 c("power", "lin", "log", "identity", "center", "hockey"))
    expect_equal(.c$covCenterType, "median")
    expect_null(.c$covCenter)
    expect_equal(.c$catCutoff, 0.05)
  })

  test_that("vaeControl validates the covariate shape settings", {
    expect_error(vaeControl(shapes = "notAShape"), "unknown covariate shape")
    ## "cat" is applied to categorical covariates automatically, never chosen
    expect_error(vaeControl(shapes = "cat"), "unknown covariate shape")
    expect_error(vaeControl(shapes = c("lin", "lin")), "duplicate")
    expect_error(vaeControl(shapes = list("power")), "named by covariate")
    expect_error(vaeControl(covCenterType = "mode"))
    ## a centering override must be a NAMED numeric
    expect_error(vaeControl(covCenter = 70))
    expect_error(vaeControl(covCenter = c(WT = NA_real_)))
    expect_error(vaeControl(catCutoff = 1.5))
    expect_error(vaeControl(catCutoff = -0.1))
  })

  test_that("the covariate shape settings round-trip through the control", {
    ## getValidNlmixrCtl.vae re-runs do.call(vaeControl, .ctl), so an argument
    ## missing from the returned list is silently dropped
    .c <- vaeControl(shapes = list(list(var = "cl", covar = "wt", shapes = "power")),
                     covCenterType = "mean", covCenter = c(WT = 70),
                     catCutoff = 0.1)
    .r <- do.call(vaeControl, .c)
    expect_equal(.r$shapes, .c$shapes)
    expect_equal(.r$covCenterType, "mean")
    expect_equal(.r$covCenter, c(WT = 70))
    expect_equal(.r$catCutoff, 0.1)
    expect_s3_class(.r, "vaeControl")
    ## and through the dispatch helper the estimation path actually uses
    .v <- getValidNlmixrCtl.vae(list(.c))
    expect_equal(.v$covCenter, c(WT = 70))
    expect_equal(.v$shapes, .c$shapes)
  })

  test_that("rxUiDeparse prints only non-default covariate shape settings", {
    .d <- rxode2::rxUiDeparse(vaeControl(), "ctl")
    expect_false(any(grepl("covCenterType", deparse(.d))))
    .n <- rxode2::rxUiDeparse(vaeControl(covCenterType = "mean"), "ctl")
    expect_true(any(grepl("covCenterType", deparse(.n))))
  })
})
