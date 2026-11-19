test_that("normal etas pass and skewed etas fail the shapiro test", {
  .x <- stats::qnorm(stats::ppoints(200))
  .r <- .etaNormalityRow("eta.cl", .x)
  expect_identical(names(.r), c("eta", "n", "mean", "sd", "meanP", "shapiroW", "shapiroP"))
  expect_identical(.r$n, 200L)
  expect_gt(.r$shapiroP, 0.5)
  expect_gt(.r$meanP, 0.5)
  # heavily skewed (log-normal) and shifted away from zero
  .y <- exp(2 * .x)
  .r <- .etaNormalityRow("eta.cl", .y)
  expect_lt(.r$shapiroP, 1e-6)
  expect_lt(.r$meanP, 1e-3)
  expect_equal(.r$shapiroW, unname(stats::shapiro.test(.y)$statistic))
  expect_equal(.r$meanP, stats::t.test(.y)$p.value)
})

test_that("small, huge, constant and skipped etas give NA statistics", {
  .na <- c("meanP", "shapiroW", "shapiroP")
  .r <- .etaNormalityRow("eta.cl", c(0.1, -0.1))
  expect_identical(.r$n, 2L)
  expect_true(all(is.na(.r[, c("mean", "sd", .na)])))
  .r <- .etaNormalityRow("eta.cl", stats::qnorm(stats::ppoints(5001)))
  expect_false(is.na(.r$meanP))
  expect_true(all(is.na(.r[, c("shapiroW", "shapiroP")])))
  .r <- .etaNormalityRow("eta.cl", rep(0, 10))
  expect_true(all(is.na(.r[, .na])))
  .r <- .etaNormalityRow("eta.cl", stats::qnorm(stats::ppoints(20)), skip = TRUE)
  expect_true(all(is.na(.r[, c("mean", "sd", .na)])))
})

test_that("only eta columns of a data frame are tested", {
  .x <- stats::qnorm(stats::ppoints(30))
  .d <- data.frame(ID = 1:30, eta.ka = .x, eta.cl = rev(.x), mixnum = 1L, occ = 1:30)
  .r <- .etaNormalityDf(.d, skip = "eta.cl", drop = "occ")
  expect_identical(.r$eta, c("eta.ka", "eta.cl"))
  expect_false(is.na(.r$shapiroP[1]))
  expect_true(is.na(.r$shapiroP[2]))
  expect_null(.etaNormalityDf(NULL))
})

test_that("fixed and zero-variance etas are found from iniDf", {
  .ui <- rxode2::rxode2(one.compartment)
  .ui <- suppressMessages(rxode2::ini(.ui, eta.ka ~ fix(0.6)))
  .ui <- suppressMessages(rxode2::ini(.ui, eta.v ~ 0))
  expect_identical(.etaNormalityFixed(.ui$iniDf), c("eta.ka", "eta.v"))
  expect_identical(.etaNormalityFixed(NULL), character(0))
})

nmTest({
  test_that("fit$etaNormality is available on a cached fit", {
    .fit <- one.compartment.fit.focei
    .r <- .fit$etaNormality
    expect_identical(.r$eta, c("eta.ka", "eta.cl", "eta.v"))
    expect_true(all(.r$n == 12L))
    expect_equal(.r$shapiroP[2], stats::shapiro.test(.fit$eta$eta.cl)$p.value)
    expect_true("etaNormality" %in% names(.nmObjGetSupportedDollars()))
  })

  test_that("a fixed-omega eta gives NA in a fit", {
    .mod <- suppressMessages(rxode2::ini(one.compartment, eta.ka ~ fix(0.6)))
    .fit <- suppressWarnings(.nlmixr(.mod, nlmixr2data::theo_sd, est = "focei", control = foceiControlFast))
    .r <- .fit$etaNormality
    expect_identical(.r$eta, c("eta.ka", "eta.cl", "eta.v"))
    expect_true(is.na(.r$shapiroP[1]))
    expect_false(is.na(.r$shapiroP[2]))
  })

  test_that("occasion etas get their own row", {
    .fit <- .iovFit("focei")
    .r <- .fit$etaNormality
    expect_identical(.r$eta, c("eta.v", "iov.cl"))
    expect_identical(.r$n, c(12L, 24L))
    expect_equal(.r$shapiroP[2], stats::shapiro.test(.fit$iov$occ$iov.cl)$p.value)
  })
})
