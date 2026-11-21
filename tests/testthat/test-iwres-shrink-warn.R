test_that(".iwresShrinkWarn() warns only for strongly negative IWRES shrinkage", {
  .mk <- function(v) {
    data.frame(
      eta.ka = c(0, 1, 1, 0, 0, 10, 10, 0, 1),
      IWRES = c(0, 1, 1, 0, 0, 2 * v, v, 0, 1),
      row.names = c(
        "mean",
        "var",
        "sd",
        "skewness",
        "kurtosis",
        "var shrinkage (%)",
        "sd shrinkage (%)",
        "t statistic",
        "p-value"
      )
    )
  }
  expect_warning(
    expect_true(.iwresShrinkWarn(.mk(-20))),
    "IWRES shrinkage is strongly negative"
  )
  expect_warning(expect_false(.iwresShrinkWarn(.mk(-4))), NA)
  expect_warning(expect_false(.iwresShrinkWarn(.mk(NA_real_))), NA)
  expect_warning(expect_false(.iwresShrinkWarn(.mk(-20)[, "eta.ka", drop = FALSE])), NA)
  expect_warning(expect_false(.iwresShrinkWarn(NULL)), NA)
  # IOV fits keep the subject-level table as `id`
  expect_warning(.iwresShrinkWarn(list(id = .mk(-20))), "IWRES shrinkage")
  expect_lt(nchar("IWRES shrinkage is strongly negative; check the residual model"), 75)
})

test_that("a fit with a too-small sigma notes IWRES shrinkage in $runInfo", {
  skip_on_cran()
  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- fix(0.05)
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  fit <- suppressMessages(nlmixr2(
    one.cmt,
    nlmixr2data::theo_sd,
    "focei",
    control = foceiControl(print = 0, maxOuterIterations = 0, covMethod = "")
  ))
  expect_lt(fit$shrink["sd shrinkage (%)", "IWRES"], .iwresShrinkWarnThreshold)
  expect_true(any(grepl("IWRES shrinkage is strongly negative", fit$runInfo)))

  one.cmt2 <- one.cmt |> ini(add.sd = 0.7)
  fit2 <- suppressMessages(nlmixr2(
    one.cmt2,
    nlmixr2data::theo_sd,
    "focei",
    control = foceiControl(print = 0, covMethod = "")
  ))
  expect_false(any(grepl("IWRES shrinkage", fit2$runInfo)))
})
