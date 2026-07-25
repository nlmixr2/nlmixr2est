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
    ## a continuous covariate contributes one column per shape family
    expect_equal(res$covariate, c("WT_power", "WT_lin", "SEX"))
    expect_equal(res$raw, c("WT", "WT", "SEX"))
    expect_equal(res$shape, c("power", "lin", "cat"))
    expect_equal(res$type, c("continuous", "continuous", "categorical"))
    ## WT is centered at its median; SEX is a 0/1 indicator so it is left RAW
    ## (center 0) -- its coefficient is the level-1 shift and the structural
    ## theta stays the reference (SEX=0) value
    expect_equal(res$center, c(rep(stats::median(c(70, 80, 60, 75)), 2), 0))
    ## the two WT shapes compete for one slot; SEX is its own group
    expect_equal(res$group, c(1L, 1L, 2L))
  })

  test_that("restricting shapes= restricts the candidate columns", {
    d <- data.frame(
      id = rep(1:4, each = 3), time = rep(0:2, 4), dv = 1:12,
      wt = rep(c(70, 80, 60, 75), each = 3))
    res <- vaeCovariates(d, shapes = "power", covCenterType = "mean")
    ## a single shape family reproduces the historic one-column-per-covariate
    ## design, mean-centered
    expect_equal(res$covariate, "WT_power")
    expect_equal(res$center, mean(c(70, 80, 60, 75)))
    ## a linear-family shape names its column after the shape asked for
    expect_equal(vaeCovariates(d, shapes = "identity")$covariate, "WT_identity")
  })

  test_that("covCenter overrides the centering statistic per covariate", {
    d <- data.frame(
      id = rep(1:4, each = 3), time = rep(0:2, 4), dv = 1:12,
      wt = rep(c(70, 80, 60, 75), each = 3))
    expect_equal(unique(vaeCovariates(d, covCenter = c(wt = 70))$center), 70)
    ## an unmatched name leaves the statistic in place
    expect_equal(unique(vaeCovariates(d, covCenter = c(age = 40))$center),
                 stats::median(c(70, 80, 60, 75)))
  })

  test_that("a non-0/1 two-level covariate becomes an indicator", {
    d <- data.frame(
      id = rep(1:4, each = 3), time = rep(0:2, 4), dv = 1:12,
      grp = rep(c(1, 2, 1, 2), each = 3))       # 2 levels, but not 0/1
    res <- vaeCovariates(d)
    expect_equal(res$type, "categorical")
    ## reference is the modal level (a tie resolves to the first alphabetically,
    ## here "1"), so the indicator is for level 2 and is never centered
    expect_equal(res$covariate, "GRP_2")
    expect_equal(res$level, "2")
    expect_equal(res$center, 0)
  })

  test_that("a factor covariate expands to indicators off the modal level", {
    d <- data.frame(
      id = rep(1:10, each = 2), time = rep(0:1, 10), dv = 1:20,
      race = rep(c("W", "W", "W", "W", "W", "B", "B", "B", "A", "A"), each = 2),
      stringsAsFactors = FALSE)
    res <- vaeCovariates(d)
    ## W holds 5/10 subjects -> reference; B and A each clear catCutoff
    expect_equal(res$covariate, c("RACE_B", "RACE_A"))
    expect_equal(res$level, c("B", "A"))
    ## every level is its own exclusion group -- several levels may enter together
    expect_equal(length(unique(res$group)), 2L)
  })

  test_that("levels below catCutoff are lumped with the reference", {
    d <- data.frame(
      id = rep(1:10, each = 2), time = rep(0:1, 10), dv = 1:20,
      race = rep(c("W", "W", "W", "W", "W", "W", "B", "B", "B", "A"), each = 2),
      stringsAsFactors = FALSE)
    ## A holds 1/10 = 0.10 of subjects, below a 0.2 cutoff
    res <- vaeCovariates(d, catCutoff = 0.2)
    expect_equal(res$covariate, "RACE_B")
    ## with no cutoff both non-reference levels are tested
    expect_equal(vaeCovariates(d, catCutoff = 0)$covariate,
                 c("RACE_B", "RACE_A"))
  })

  test_that("log shapes are skipped for a non-positive covariate", {
    d <- data.frame(
      id = rep(1:4, each = 3), time = rep(0:2, 4), dv = 1:12,
      score = rep(c(-1, 0, 2, 3), each = 3))
    res <- vaeCovariates(d)
    ## log(score/ctr) is undefined, so only the linear family survives
    expect_false(any(res$shape %in% c("power", "log")))
    expect_true(all(res$type == "continuous"))
    ## and asking ONLY for a log shape still leaves the covariate searchable
    expect_equal(nrow(vaeCovariates(d, shapes = "power")), 1L)
  })

  test_that("vaeCovariates matches what .vaeDataPrep explores (theo_sd)", {
    res <- vaeCovariates(nlmixr2data::theo_sd)
    expect_equal(res$raw, c("WT", "WT"))
    expect_equal(res$type, c("continuous", "continuous"))
  })

  test_that("vaeCovariates warns on time-varying columns and drops them", {
    d <- data.frame(
      id = rep(1:2, each = 3),
      time = rep(0:2, 2),
      dv = 1:6,
      crcl = c(90, 85, 80, 100, 95, 90),   # varies within subject
      wt = rep(c(70, 80), each = 3))
    expect_warning(res <- vaeCovariates(d), "time-varying")
    expect_equal(unique(res$raw), "WT")
    ## warn=FALSE drops them silently
    expect_silent(res2 <- vaeCovariates(d, warn = FALSE))
    expect_equal(res2, res)
  })

  test_that("vaeCovariates returns zero rows when nothing qualifies", {
    d <- data.frame(id = rep(1:2, each = 2), time = rep(0:1, 2), dv = 1:4)
    res <- vaeCovariates(d)
    expect_equal(nrow(res), 0L)
    expect_equal(names(res), c("covariate", "raw", "shape", "level", "group",
                               "type", "center"))
  })

  test_that("vaeCovariates requires an ID column", {
    expect_error(vaeCovariates(data.frame(time = 0:1, dv = 1:2)), "ID")
  })
})
