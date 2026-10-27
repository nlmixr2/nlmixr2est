# Every estimator that takes IOV fits a model whose only random effect is an
# occasion one, in each IOV expansion
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
    linCmt() ~ add(add.sd)
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
      list("focei", "theta", 200.866286),
      list("focei", "omega", 200.312479),
      list("foce", "theta", 201.415327),
      list("foce", "omega", 200.294408)
    )
    for (.w in .want) {
      .fit <- suppressWarnings(suppressMessages(nlmixr2(
        .iovOnlyFitMod,
        .iovOnlyFitData(),
        .w[[1]],
        control = list(print = 0, covMethod = "", iovMethod = .w[[2]])
      )))
      .label <- paste(.w[[1]], .w[[2]])
      expect_equal(.fit$objf, .w[[3]], tolerance = 1e-6, label = .label)
      .iovOnlyFitShape(.fit, .label)
    }
    # saem fits it; its objective is finite (which objective it reports is
    # its quadrature's question, not this one's)
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.iovOnlyFitMod, .iovOnlyFitData(), "saem", control = list(print = 0, covMethod = ""))
    ))
    expect_true(is.finite(suppressMessages(.fit$objf)))
    .iovOnlyFitShape(.fit, "saem")
  })

  test_that("laplace, agq, imp and npag fit a model with only IOV", {
    # laplace and agq are deterministic, as focei is; imp samples and npag
    # places support points, so they are held to the fit's shape
    .want <- c(laplace = 200.866286, agq = 198.7465, imp = NA, npag = NA)
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
        expect_equal(.fit$objf, .want[[.est]], tolerance = 1e-6, label = .est)
      }
      .iovOnlyFitShape(.fit, .est)
    }
  })
})
