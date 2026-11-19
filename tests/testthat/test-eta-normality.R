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
  # positive sd, but t.test() calls it essentially constant
  .r <- .etaNormalityRow("eta.cl", c(1, 1, 1 + .Machine$double.eps))
  expect_true(is.na(.r$meanP))
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

test_that("component-owned mixture etas use only their component's subjects", {
  .ui <- rxode2::rxode2(function() {
    ini({
      tcl1 <- 0
      tcl2 <- 2
      p1 <- 0.5
      tv <- 3
      eta.cl1 ~ 0.1
      eta.cl2 ~ 0.1
      eta.v ~ 0.1
      add.sd <- 0.1
    })
    model({
      cl <- mix(exp(tcl1 + eta.cl1), p1, exp(tcl2 + eta.cl2))
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  })
  .x <- stats::qnorm(stats::ppoints(10))
  .mn <- rep(1:2, each = 10)
  .ranef <- data.frame(
    ID = 1:20,
    eta.cl1 = c(.x, rep(0, 10)),
    eta.cl2 = c(rep(0, 10), .x),
    eta.v = c(.x, .x),
    mixnum = .mn
  )
  .r <- .etaNormalityFit(list(ui = .ui, ranef = .ranef, mixNum = data.frame(ID = 1:20, mixnum = .mn)))
  expect_identical(.r$eta, c("eta.cl1", "eta.cl2", "eta.v"))
  expect_identical(.r$n, c(10L, 10L, 20L))
  expect_equal(.r$shapiroP[1], stats::shapiro.test(.x)$p.value)
  # without the component, the placeholder zeros are counted
  .r <- .etaNormalityFit(list(ui = .ui, ranef = .ranef))
  expect_identical(.r$n, c(20L, 20L, 20L))
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

  test_that("a zero-omega eta dropped by the fit still gets an NA row", {
    .mod <- suppressMessages(rxode2::ini(one.compartment, eta.ka ~ 0))
    .fit <- suppressWarnings(.nlmixr(.mod, nlmixr2data::theo_sd, est = "focei", control = foceiControlFast))
    .r <- .fit$etaNormality
    expect_setequal(.r$eta, c("eta.ka", "eta.cl", "eta.v"))
    .ka <- .r[.r$eta == "eta.ka", ]
    expect_true(all(is.na(.ka[, c("mean", "sd", "meanP", "shapiroW", "shapiroP")])))
    expect_false(anyNA(.r$shapiroP[.r$eta != "eta.ka"]))
  })

  test_that("occasion etas get their own row", {
    .fit <- .iovFit("focei")
    .r <- .fit$etaNormality
    expect_identical(.r$eta, c("eta.v", "iov.cl"))
    expect_identical(.r$n, c(12L, 24L))
    expect_equal(.r$shapiroP[2], stats::shapiro.test(.fit$iov$occ$iov.cl)$p.value)
  })
})
