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
  expect_true(all(is.na(.etaNormalityShapiro(NULL))))
  expect_true(all(is.na(.etaNormalityShapiro(c(0.1, -0.1)))))
  expect_true(all(is.na(.etaNormalityShapiro(stats::qnorm(stats::ppoints(5001))))))
  expect_true(all(is.na(.etaNormalityShapiro(rep(0, 10)))))
})

test_that("shapiro rows are appended to the shrinkage table", {
  .x <- stats::qnorm(stats::ppoints(20))
  .s <- .etaNormalityShrinkTab(eta.cl = c(0.1, 0.3, 0.04), IWRES = c(0, 1, 1))
  .r <- .shrinkAddShapiro(.s, list(eta.cl = .x))
  expect_identical(row.names(.r), c(row.names(.s), .shrinkShapiroRows))
  expect_equal(.r[1:9, ], .s)
  expect_equal(.r["Shapiro-Wilk p-value", "eta.cl"], stats::shapiro.test(.x)$p.value)
  expect_true(is.na(.r["Shapiro-Wilk W", "IWRES"]))
  # adding twice keeps one set of rows
  expect_identical(.shrinkAddShapiro(.r, list(eta.cl = .x)), .r)
  # IOV fits keep one table per level
  .l <- list(id = .s, occ = .etaNormalityShrinkTab(iov.cl = c(-0.2, 0.5, 0.6)))
  .r <- .shrinkAddShapiro(.l, list(eta.cl = .x, iov.cl = rev(.x)))
  expect_identical(names(.r), c("id", "occ"))
  expect_equal(.r$occ["Shapiro-Wilk W", "iov.cl"], unname(stats::shapiro.test(.x)$statistic))
  expect_null(.shrinkAddShapiro(NULL, list()))
})

test_that("an eta's statistics are read from the shrinkage table", {
  .x <- stats::qnorm(stats::ppoints(20))
  .s <- .shrinkAddShapiro(.etaNormalityShrinkTab(eta.cl = c(0.1, 0.3, 0.04)), list(eta.cl = .x))
  .r <- .etaNormalityShrink(.s, "eta.cl")
  expect_identical(names(.r), c("mean", "sd", "meanP", "shapiroW", "shapiroP"))
  expect_equal(.r[1:3], c(mean = 0.1, sd = 0.3, meanP = 0.04))
  expect_equal(.r[["shapiroP"]], stats::shapiro.test(.x)$p.value)
  expect_true(all(is.na(.etaNormalityShrink(.s, "eta.v"))))
  expect_true(all(is.na(.etaNormalityShrink(NULL, "eta.cl"))))
  .l <- list(id = .s, occ = .etaNormalityShrinkTab(iov.cl = c(-0.2, 0.5, 0.6)))
  expect_equal(.etaNormalityShrink(.l, "iov.cl")[1:3], c(mean = -0.2, sd = 0.5, meanP = 0.6))
})

test_that("fixed, zero-variance and dropped etas give NA rows", {
  .x <- stats::qnorm(stats::ppoints(20))
  .values <- list(eta = list(eta.ka = .x, eta.cl = .x), occ = list())
  .s <- .shrinkAddShapiro(.etaNormalityShrinkTab(eta.ka = c(0, 1, 1), eta.cl = c(0, 1, 1)), .values$eta)
  .r <- .etaNormalityDf(.values, .s, skip = c("eta.cl", "eta.v"))
  expect_identical(names(.r), c("eta", "n", "mean", "sd", "meanP", "shapiroW", "shapiroP"))
  expect_identical(.r$eta, c("eta.ka", "eta.cl", "eta.v"))
  expect_identical(.r$n, c(20L, 20L, NA_integer_))
  expect_false(anyNA(.r[1, -1]))
  expect_true(all(is.na(.r[2:3, 3:7])))
  expect_null(.etaNormalityDf(list(eta = list()), NULL))
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
  .mixNum <- data.frame(ID = 1:20, mixnum = .mn)
  .v <- .etaNormalityValues(list(ui = .ui, ranef = .ranef, mixNum = .mixNum))
  expect_identical(names(.v$eta), c("eta.cl1", "eta.cl2", "eta.v"))
  expect_equal(.v$eta$eta.cl1, .x)
  expect_equal(.v$eta$eta.cl2, .x)
  expect_equal(.v$eta$eta.v, c(.x, .x))
  # the component can also come from ranef's own mixnum column
  expect_equal(.etaNormalityValues(list(ui = .ui, ranef = .ranef)), .v)
  # without any component, the placeholder zeros are counted
  .v <- .etaNormalityValues(list(ui = .ui, ranef = .ranef[, names(.ranef) != "mixnum"]))
  expect_equal(lengths(.v$eta), c(eta.cl1 = 20L, eta.cl2 = 20L, eta.v = 20L))
  # an owned occasion eta is filtered by its subject's component, and split
  # by occasion for the per-occasion shrinkage columns
  .own <- c(.x, -.x)
  .iov <- data.frame(ID = rep(1:20, each = 2), occ = rep(1:2, 20), eta.cl2 = c(rep(0, 20), .own))
  .v <- .etaNormalityValues(list(
    ui = .ui,
    ranef = .ranef[, c("ID", "eta.v", "mixnum")],
    iov = list(occ = .iov)
  ))
  expect_identical(names(.v$eta), c("eta.v", "eta.cl2"))
  expect_equal(.v$eta$eta.cl2, .own)
  expect_identical(names(.v$occ), c("rx.eta.cl2.1", "rx.eta.cl2.2"))
  expect_equal(.v$occ$rx.eta.cl2.1, .own[c(TRUE, FALSE)])
})

nmTest({
  test_that("$shrink and $etaNormality agree on a cached fit", {
    .fit <- one.compartment.fit.focei
    .s <- .fit$shrink
    expect_identical(tail(row.names(.s), 2), .shrinkShapiroRows)
    expect_equal(.s["Shapiro-Wilk p-value", "eta.cl"], stats::shapiro.test(.fit$eta$eta.cl)$p.value)
    expect_true(is.na(.s["Shapiro-Wilk p-value", "IWRES"]))
    # the stored table is left as it was
    expect_false(any(.shrinkShapiroRows %in% row.names(get("shrink", envir = .fit$env, inherits = FALSE))))
    .r <- .fit$etaNormality
    expect_identical(.r$eta, c("eta.ka", "eta.cl", "eta.v"))
    expect_true(all(.r$n == 12L))
    expect_equal(.r$mean, unname(unlist(.s["mean", .r$eta])))
    expect_equal(.r$sd, unname(unlist(.s["sd", .r$eta])))
    expect_equal(.r$meanP, unname(unlist(.s["p-value", .r$eta])))
    expect_equal(.r$shapiroW, unname(unlist(.s["Shapiro-Wilk W", .r$eta])))
    expect_equal(.r$shapiroP, unname(unlist(.s["Shapiro-Wilk p-value", .r$eta])))
    expect_true(all(c("etaNormality", "shrink") %in% names(.nmObjGetSupportedDollars())))
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

  test_that("occasion etas get shapiro rows per occasion and pooled", {
    .fit <- .iovFit("focei")
    .s <- .fit$shrink
    .iov <- .fit$iov$occ
    expect_equal(.s$occ["Shapiro-Wilk p-value", "iov.cl"], stats::shapiro.test(.iov$iov.cl)$p.value)
    expect_equal(
      .s$occ["Shapiro-Wilk p-value", "rx.iov.cl.1"],
      stats::shapiro.test(.iov$iov.cl[.iov$occ == 1])$p.value
    )
    .r <- .fit$etaNormality
    expect_identical(.r$eta, c("eta.v", "iov.cl"))
    expect_identical(.r$n, c(12L, 24L))
    expect_equal(.r$meanP, c(.s$id["p-value", "eta.v"], .s$occ["p-value", "iov.cl"]))
    expect_equal(.r$shapiroP[2], .s$occ["Shapiro-Wilk p-value", "iov.cl"])
  })
})
