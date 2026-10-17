# A model whose only random effect is an occasion one: no subject-level eta
# for the IOV rewrite to copy its rows from
.iovOnlyMod <- function() {
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

.iovOnlyData <- function() {
  .d <- nlmixr2data::theo_sd
  .d$occ <- ifelse(.d$TIME < 6, 1L, 2L)
  .d
}

nmTest({
  test_that("the IOV rewrite of a model with no subject-level eta keeps whole rows", {
    .ui <- rxode2::rxode2(.iovOnlyMod)
    # theta expansion: a magnitude theta and unit occasion etas at "id"
    .idf <- suppressMessages(.uiApplyIov(.ui, "focei", .iovOnlyData(), foceiControl(iovMethod = "theta")))$ui$iniDf
    .th <- .idf[.idf$name == "iov.cl", ]
    expect_identical(.th$ntheta, 5L)
    expect_identical(.th$lower, 0)
    expect_equal(.th$est, sqrt(0.1))
    expect_false(.th$fix)
    .eta <- .idf[!is.na(.idf$neta1), ]
    expect_identical(.eta$name, c("rx.iov.cl.1", "rx.iov.cl.2"))
    expect_identical(.eta$condition, c("id", "id"))
    expect_identical(.eta$est, c(1, 1))
    expect_identical(.eta$fix, c(TRUE, TRUE))
    # omega expansion: the occasion etas carry the variance, the second a
    # same() copy of the first
    .idf <- suppressMessages(.uiApplyIov(.ui, "focei", .iovOnlyData(), foceiControl(iovMethod = "omega")))$ui$iniDf
    .eta <- .idf[!is.na(.idf$neta1), ]
    expect_identical(.eta$condition, c("id", "id:same:rx.iov.cl.1"))
    expect_identical(.eta$est, c(0.1, 0.1))
  })

  test_that("a focei fit with only IOV keeps its eta and shrinkage tables", {
    skip_on_cran()
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.iovOnlyMod, .iovOnlyData(), "focei", foceiControl(print = 0, covMethod = ""))
    ))
    expect_equal(.fit$objf, 200.866286, tolerance = 1e-6)
    # no subject-level eta is left once the occasion columns move out
    expect_s3_class(.fit$eta, "data.frame")
    expect_identical(names(.fit$eta), "ID")
    expect_identical(nrow(.fit$eta), 12L)
    expect_identical(typeof(.fit$etaMat), "double")
    expect_identical(dim(.fit$etaMat), c(12L, 2L))
    expect_s3_class(.fit$shrink$id, "data.frame")
  })
})
