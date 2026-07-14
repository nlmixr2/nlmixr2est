## vaeCovariates(): exported view of the VAE automatic covariate search
## (fast, data-only -- no model build or training)

nmTest({
  test_that("vaeCovariates discovers subject-level covariates with fit rules", {
    d <- data.frame(
      id = rep(1:4, each = 3),
      time = rep(0:2, 4),
      dv = 1:12,
      wt = rep(c(70, 80, 60, 75), each = 3),
      sex = rep(c(0, 1, 0, 1), each = 3))
    res <- vaeCovariates(d)
    expect_s3_class(res, "data.frame")
    expect_equal(res$covariate, c("WT", "SEX"))
    expect_equal(res$type, c("continuous", "categorical"))
    expect_equal(res$center, c(mean(c(70, 80, 60, 75)), 0.5))
  })

  test_that("vaeCovariates matches what .vaeDataPrep explores (theo_sd)", {
    res <- vaeCovariates(nlmixr2data::theo_sd)
    expect_equal(res$covariate, "WT")
    expect_equal(res$type, "continuous")
  })

  test_that("vaeCovariates warns on time-varying columns and drops them", {
    d <- data.frame(
      id = rep(1:2, each = 3),
      time = rep(0:2, 2),
      dv = 1:6,
      crcl = c(90, 85, 80, 100, 95, 90),   # varies within subject
      wt = rep(c(70, 80), each = 3))
    expect_warning(res <- vaeCovariates(d), "time-varying")
    expect_equal(res$covariate, "WT")
    ## warn=FALSE drops them silently
    expect_silent(res2 <- vaeCovariates(d, warn = FALSE))
    expect_equal(res2, res)
  })

  test_that("vaeCovariates returns zero rows when nothing qualifies", {
    d <- data.frame(id = rep(1:2, each = 2), time = rep(0:1, 2), dv = 1:4)
    res <- vaeCovariates(d)
    expect_equal(nrow(res), 0L)
    expect_equal(names(res), c("covariate", "type", "center"))
  })

  test_that("vaeCovariates requires an ID column", {
    expect_error(vaeCovariates(data.frame(time = 0:1, dv = 1:2)), "ID")
  })
})
