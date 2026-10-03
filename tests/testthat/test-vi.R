# Installing the population variational covariance of a full-Bayes vi fit
# (R/vi.R), without a fit.

.viRes <- function(viCov) {
  list(phiThetaIdx = c(0L, 1L, -1L), prep = list(thetaRealNames = c("tka", "add.sd")), viCov = viCov)
}

.viEnv <- function() {
  .pf <- data.frame(
    Parameter = c("", ""),
    Estimate = c(0.45, 0.7),
    SE = c(NA_real_, NA_real_),
    "%RSE" = c(NA_real_, NA_real_),
    "Back-transformed" = c(exp(0.45), 0.7),
    "CI Lower" = c(NA_real_, NA_real_),
    "CI Upper" = c(NA_real_, NA_real_),
    check.names = FALSE,
    row.names = c("tka", "add.sd")
  )
  .e <- new.env(parent = emptyenv())
  .e$parFixedDf <- .pf
  .e$parFixed <- .updateParFixedApplySig(.pf, 3L, 0.9, character())
  class(.e$parFixed) <- c("nlmixr2ParFixed", "data.frame")
  .e$control <- list(ci = 0.9, sigdigTable = 3L)
  .e$covMethod <- ""
  .e
}

.viCov <- matrix(c(0.04, 0.002, 0.5, 0.002, 0.01, 0.1, 0.5, 0.1, 9), 3)

test_that("the variational covariance refreshes the formatted table, with the identity interval", {
  .e <- .viEnv()
  .adviInstallVarCov(.e, .viRes(.viCov))
  expect_identical(.e$covMethod, "vi")
  expect_equal(unname(.e$cov), .viCov[1:2, 1:2])
  expect_identical(rownames(.e$cov), c("tka", "add.sd"))
  expect_equal(.e$parFixedDf$SE, c(0.2, 0.1))
  expect_equal(.e$parFixedDf[["%RSE"]], c(0.2 / 0.45, 0.1 / 0.7) * 100)
  # the identity interval, only where the estimate is reported untransformed
  .z <- stats::qnorm(0.95)
  expect_equal(.e$parFixedDf["add.sd", "CI Lower"], 0.7 - .z * 0.1)
  expect_equal(.e$parFixedDf["add.sd", "CI Upper"], 0.7 + .z * 0.1)
  expect_true(is.na(.e$parFixedDf["tka", "CI Lower"]))
  # the printed table follows
  expect_equal(as.numeric(.e$parFixed[, "SE"]), c(0.2, 0.1))
  expect_true("Back-transformed(90%CI)" %in% names(.e$parFixed))
})

test_that("a variational covariance that is not finite is not installed", {
  .e <- .viEnv()
  .bad <- .viCov
  .bad[2, 2] <- NaN
  expect_warning(
    .adviInstallVarCov(.e, .viRes(.bad)),
    "\"vi\" covariance is not finite; none installed",
    fixed = TRUE
  )
  expect_identical(.e$covMethod, "")
  expect_false(exists("cov", envir = .e, inherits = FALSE))
  expect_true(all(is.na(.e$parFixedDf$SE)))
})
