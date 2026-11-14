# Every estimator that takes IOV fits a model whose only random effect is an
# occasion one, in each IOV expansion.  The model is an ODE model: occ steps
# within a subject, and how linCmt() reads a covariate over an interval
# depends on the rxode2 version.  The objectives move by ~1e-5 relative
# between platforms (CI vs a local build), so they are pinned to 1e-4.
.iovOnlyFitMod <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    iov.cl ~ 0.1 | occ
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + iov.cl)
    v <- exp(tv)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

.iovOnlyFitData <- function() {
  .d <- nlmixr2data::theo_sd
  .d$occ <- ifelse(.d$TIME < 6, 1L, 2L)
  .d
}

# what an IOV-only fit keeps: no subject-level eta, and the occasion
# estimates by subject and occasion
.iovOnlyFitShape <- function(fit, label) {
  expect_identical(names(fit$eta), "ID", label = label)
  expect_identical(typeof(fit$etaMat), "double", label = label)
  expect_identical(names(fit$ranef), "ID", label = label)
  expect_identical(nrow(fit$ranef), 12L, label = label)
  .iov <- fit$iov$occ
  expect_identical(names(.iov), c("ID", "occ", "iov.cl"), label = label)
  expect_identical(nrow(.iov), 24L, label = label)
  expect_true(all(is.finite(.iov$iov.cl)), label = label)
}

nmTest({
  test_that("focei, foce and saem fit a model with only IOV", {
    .want <- list(
      list("focei", "theta", 193.607959),
      list("focei", "omega", 193.608374),
      list("foce", "theta", 193.608814),
      list("foce", "omega", 193.609741)
    )
    for (.w in .want) {
      .fit <- suppressWarnings(suppressMessages(nlmixr2(
        .iovOnlyFitMod,
        .iovOnlyFitData(),
        .w[[1]],
        control = list(print = 0, covMethod = "", iovMethod = .w[[2]])
      )))
      .label <- paste(.w[[1]], .w[[2]])
      expect_equal(.fit$objf, .w[[3]], tolerance = 1e-4, label = .label)
      .iovOnlyFitShape(.fit, .label)
    }
    # saem fits it, and reports its Gaussian-quadrature objective
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.iovOnlyFitMod, .iovOnlyFitData(), "saem", control = list(print = 0, covMethod = ""))
    ))
    expect_equal(suppressMessages(.fit$objf), 197.935147, tolerance = 1e-4)
    expect_identical(.fit$ofvType, "gauss3_1.6")
    .iovOnlyFitShape(.fit, "saem")
  })

  test_that("saem's two-level IOV methods fall back to the theta expansion", {
    # two-level and collapsed IOV ride on a subject-level eta of the same
    # parameter, which an IOV-only model does not have
    for (.m in c("twoLevel", "collapsed", "theta")) {
      .fit <- suppressWarnings(suppressMessages(nlmixr2(
        .iovOnlyFitMod,
        .iovOnlyFitData(),
        "saem",
        control = list(print = 0, covMethod = "", iovMethod = .m)
      )))
      expect_equal(suppressMessages(.fit$objf), 197.935147, tolerance = 1e-4, label = .m)
      expect_identical(
        any(grepl("used iovMethod='theta'", .fit$runInfo, fixed = TRUE)),
        .m != "theta",
        label = .m
      )
      # the magnitude's bound is the rewrite's own, so saem has none to report
      expect_false(any(grepl("user-defined boundaries", .fit$runInfo, fixed = TRUE)), label = .m)
      .iovOnlyFitShape(.fit, paste("saem", .m))
    }
  })

  test_that("focei fits a model with only IOV under every iovXform", {
    # var stops ~0.2 above the others, at a point that moves by platform
    .want <- c(sd = 193.607959, var = 193.832538, logsd = 193.608392, logvar = 193.608758)
    .tol <- c(sd = 1e-4, var = 2e-3, logsd = 1e-4, logvar = 1e-4)
    for (.x in names(.want)) {
      .fit <- suppressWarnings(suppressMessages(nlmixr2(
        .iovOnlyFitMod,
        .iovOnlyFitData(),
        "focei",
        control = list(print = 0, covMethod = "", iovMethod = "theta", iovXform = .x)
      )))
      expect_equal(.fit$objf, .want[[.x]], tolerance = .tol[[.x]], label = .x)
      .iovOnlyFitShape(.fit, paste("focei", .x))
    }
  })

  test_that("laplace, agq, imp and npag fit a model with only IOV", {
    # laplace and agq are deterministic, as focei is; imp samples and npag
    # places support points, so they are held to the fit's shape
    .want <- c(laplace = 193.607959, agq = 194.033970, imp = NA, npag = NA)
    for (.est in names(.want)) {
      .fit <- suppressWarnings(suppressMessages(nlmixr2(
        .iovOnlyFitMod,
        .iovOnlyFitData(),
        .est,
        control = list(print = 0, covMethod = "")
      )))
      if (is.na(.want[[.est]])) {
        expect_true(is.finite(.fit$objf), label = .est)
      } else {
        expect_equal(.fit$objf, .want[[.est]], tolerance = 1e-4, label = .est)
      }
      .iovOnlyFitShape(.fit, .est)
    }
  })
})
