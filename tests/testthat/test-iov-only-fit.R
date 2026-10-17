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
      expect_identical(names(.fit$eta), "ID", label = .label)
      expect_identical(typeof(.fit$etaMat), "double", label = .label)
    }
    # saem's objective is its Gaussian quadrature, not a FOCEi fallback
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.iovOnlyFitMod, .iovOnlyFitData(), "saem", control = list(print = 0, covMethod = ""))
    ))
    .msg <- capture.output(.objf <- .fit$objf, type = "message")
    expect_false(any(grepl("quadrature failed", .msg)))
    expect_identical(.fit$ofvType, "gauss3_1.6")
    expect_equal(.objf, 200.859435, tolerance = 1e-6)
    expect_identical(typeof(.fit$etaMat), "double")
  })
})
