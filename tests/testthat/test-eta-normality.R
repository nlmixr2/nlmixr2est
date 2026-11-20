.etaNormalityShrinkTab <- function(...) {
  .v <- list(...)
  .d <- as.data.frame(lapply(.v, function(v) c(v[1], v[2]^2, v[2], 0, 0, 0, 0, 0, v[3])))
  names(.d) <- names(.v)
  row.names(.d) <- c(
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
  .d
}

test_that("normal etas pass and skewed etas fail the shapiro test", {
  .x <- stats::qnorm(stats::ppoints(200))
  .s <- .etaNormalityShapiro(.x)
  expect_identical(names(.s), c("shapiroW", "shapiroP"))
  expect_gt(.s[["shapiroP"]], 0.5)
  # heavily skewed (log-normal)
  .y <- exp(2 * .x)
  .s <- .etaNormalityShapiro(.y)
  expect_lt(.s[["shapiroP"]], 1e-6)
  expect_equal(.s[["shapiroW"]], unname(stats::shapiro.test(.y)$statistic))
})

test_that("small, huge and constant etas have no shapiro test", {
  expect_true(all(is.na(.etaNormalityShapiro(c(0.1, -0.1)))))
  expect_true(all(is.na(.etaNormalityShapiro(stats::qnorm(stats::ppoints(5001))))))
  expect_true(all(is.na(.etaNormalityShapiro(rep(0, 10)))))
})

test_that("mean, sd and meanP come from $shrink", {
  .s <- .etaNormalityShrinkTab(eta.cl = c(0.1, 0.3, 0.04))
  expect_equal(.etaNormalityShrink(.s, "eta.cl"), c(mean = 0.1, sd = 0.3, meanP = 0.04))
  expect_true(all(is.na(.etaNormalityShrink(.s, "eta.v"))))
  expect_true(all(is.na(.etaNormalityShrink(NULL, "eta.cl"))))
  # IOV fits keep one table per level
  .l <- list(id = .s, occ = .etaNormalityShrinkTab(iov.cl = c(-0.2, 0.5, 0.6)))
  expect_equal(.etaNormalityShrink(.l, "iov.cl"), c(mean = -0.2, sd = 0.5, meanP = 0.6))
  .r <- .etaNormalityRow("eta.cl", stats::qnorm(stats::ppoints(20)), shrink = .s)
  expect_identical(names(.r), c("eta", "n", "mean", "sd", "meanP", "shapiroW", "shapiroP"))
  expect_identical(.r$n, 20L)
  expect_equal(.r$meanP, 0.04)
  expect_false(is.na(.r$shapiroP))
  .r <- .etaNormalityRow("eta.cl", stats::qnorm(stats::ppoints(20)), shrink = .s, skip = TRUE)
  expect_true(all(is.na(.r[, c("mean", "sd", "meanP", "shapiroW", "shapiroP")])))
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
  # the component can also come from ranef's own mixnum column
  .r <- .etaNormalityFit(list(ui = .ui, ranef = .ranef))
  expect_identical(.r$n, c(10L, 10L, 20L))
  # an owned occasion eta is filtered by its subject's component too
  .iov <- data.frame(ID = rep(1:20, each = 2), occ = rep(1:2, 20), eta.cl2 = c(rep(0, 20), .x, .x))
  .r <- .etaNormalityFit(list(
    ui = .ui,
    ranef = .ranef[, c("ID", "eta.v")],
    iov = list(occ = .iov),
    mixNum = data.frame(ID = 1:20, mixnum = .mn)
  ))
  expect_identical(.r$eta, c("eta.v", "eta.cl2"))
  expect_identical(.r$n, c(20L, 20L))
  expect_equal(.r$shapiroP[2], stats::shapiro.test(c(.x, .x))$p.value)
  # ... including when the component is only on ranef
  .r <- .etaNormalityFit(list(
    ui = .ui,
    ranef = .ranef[, c("ID", "eta.v", "mixnum")],
    iov = list(occ = .iov)
  ))
  expect_identical(.r$n, c(20L, 20L))
  # without any component, the placeholder zeros are counted
  .r <- .etaNormalityFit(list(ui = .ui, ranef = .ranef[, names(.ranef) != "mixnum"]))
  expect_identical(.r$n, c(20L, 20L, 20L))
})

nmTest({
  test_that("fit$etaNormality is available on a cached fit", {
    .fit <- one.compartment.fit.focei
    .r <- .fit$etaNormality
    expect_identical(.r$eta, c("eta.ka", "eta.cl", "eta.v"))
    expect_true(all(.r$n == 12L))
    expect_equal(.r$shapiroP[2], stats::shapiro.test(.fit$eta$eta.cl)$p.value)
    .s <- .fit$shrink
    expect_equal(.r$mean, unname(unlist(.s["mean", .r$eta])))
    expect_equal(.r$sd, unname(unlist(.s["sd", .r$eta])))
    expect_equal(.r$meanP, unname(unlist(.s["p-value", .r$eta])))
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
    expect_equal(.r$meanP, c(.fit$shrink$id["p-value", "eta.v"], .fit$shrink$occ["p-value", "iov.cl"]))
  })
})
