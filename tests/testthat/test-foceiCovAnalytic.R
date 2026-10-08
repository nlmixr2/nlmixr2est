# Installing the analytic covariance (R/foceiCovAnalytic.R): the in-fit
# installer and foceiCovAnalytic(), without a fit.

.anNames <- c("tka", "tcl", "om.eta.cl")
.anPd <- matrix(c(0.04, 0.01, 0.002, 0.01, 0.09, 0.003, 0.002, 0.003, 0.05), 3, dimnames = list(.anNames, .anNames))
# indefinite, while its theta block is positive definite on its own (#1055)
.anBad <- matrix(c(0.04, 0.01, 0.05, 0.01, 0.09, 0.02, 0.05, 0.02, -0.01), 3, dimnames = list(.anNames, .anNames))

.anEnv <- function(covF, covMethod = "r", covFull = TRUE) {
  .e <- new.env(parent = emptyenv())
  .e$ui <- new.env(parent = emptyenv())
  .e$ui$control <- list(covMethod = 2L, covFull = covFull)
  .e$.analyticCov <- covF
  .e$.analyticThetaNames <- c("tka", "tcl")
  # the native step inverts the theta block of the analytic information
  .e$cov <- covF[c("tka", "tcl"), c("tka", "tcl")]
  .e$covMethod <- covMethod
  .e$objDf <- data.frame(OBJF = 1, "Condition#(Cov)" = 99, "Condition#(Cor)" = 99, check.names = FALSE)
  .e
}

test_that("a rejected analytic covariance names the theta block the native step kept", {
  .e <- .anEnv(.anBad)
  expect_warning(
    .foceiInstallAnalyticCov(.e),
    "\"analytic (full)\" covariance is not positive definite; kept \"r (analytic)\"",
    fixed = TRUE
  )
  expect_identical(.e$covMethod, "r (analytic)")
  expect_equal(.e$cov, .anBad[1:2, 1:2])
  expect_false(exists("covList", envir = .e, inherits = FALSE))
  # a corrected native R keeps its decoration; an S-matrix fallback is a
  # finite-difference covariance and keeps its own name
  .e <- .anEnv(.anBad, "|r|", covFull = FALSE)
  expect_warning(
    .foceiInstallAnalyticCov(.e),
    "\"analytic\" covariance is not positive definite; kept \"|r| (analytic)\"",
    fixed = TRUE
  )
  expect_identical(.e$covMethod, "|r| (analytic)")
  .e <- .anEnv(.anBad, "s")
  expect_warning(.foceiInstallAnalyticCov(.e), "kept \"s\"", fixed = TRUE)
  expect_identical(.e$covMethod, "s")
})

test_that("a positive-definite analytic covariance installs one shape and caches the other", {
  .e <- .anEnv(.anPd)
  .foceiInstallAnalyticCov(.e)
  expect_identical(.e$covMethod, "analytic (full)")
  expect_equal(.e$cov, .anPd)
  expect_identical(names(.e$covList), "analytic")
  expect_equal(.e$covList$analytic, .anPd[1:2, 1:2])
  .ev <- eigen(.anPd, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(.e$objDf[["Condition#(Cov)"]], max(.ev) / min(.ev))
  .e <- .anEnv(.anPd, covFull = FALSE)
  .foceiInstallAnalyticCov(.e)
  expect_identical(.e$covMethod, "analytic")
  expect_equal(.e$cov, .anPd[1:2, 1:2])
  expect_equal(.e$covList[["analytic (full)"]], .anPd)
  .ev <- eigen(.anPd[1:2, 1:2], symmetric = TRUE, only.values = TRUE)$values
  expect_equal(.e$objDf[["Condition#(Cov)"]], max(.ev) / min(.ev))
})

test_that("foceiCovAnalytic() refreshes the SEs and condition numbers of what it installs", {
  .e <- .anEnv(.anPd)
  .e$covMethod <- "r,s (full)"
  .e$parFixedDf <- data.frame(
    Estimate = c(0.45, 1),
    SE = c(NA_real_, NA_real_),
    "%RSE" = c(NA_real_, NA_real_),
    check.names = FALSE,
    row.names = c("tka", "tcl")
  )
  local_mocked_bindings(.foceiCovAnalyticCalc = function(fit) {
    list(cov = .anPd, se = sqrt(diag(.anPd)), R = solve(.anPd), params = .anNames, method = "analytic", pd = TRUE)
  })
  .r <- foceiCovAnalytic(.e)
  expect_equal(.r$cov, .anPd)
  expect_identical(.e$covMethod, "analytic (full)")
  expect_equal(.e$parFixedDf$SE, c(0.2, 0.3))
  expect_equal(.e$parFixedDf[["%RSE"]], c(0.2 / 0.45, 0.3) * 100)
  .ev <- eigen(.anPd, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(.e$objDf[["Condition#(Cov)"]], max(.ev) / min(.ev))
  # the covariance it replaced stays recoverable
  expect_equal(.e$covList[["r,s (full)"]], .anPd[1:2, 1:2])
  expect_equal(.e$covList$analytic, .anPd[1:2, 1:2])
})

test_that("foceiCovAnalytic() never installs an indefinite covariance, and says so on every call", {
  .e <- .anEnv(.anBad)
  .e$covMethod <- "r,s (full)"
  local_mocked_bindings(.foceiCovAnalyticCalc = function(fit) {
    list(
      cov = .anBad,
      se = sqrt(abs(diag(.anBad))),
      R = solve(.anBad),
      params = .anNames,
      method = "analytic",
      pd = FALSE
    )
  })
  .msg <- "\"analytic (full)\" covariance is not positive definite; kept \"r,s (full)\""
  expect_warning(.r <- foceiCovAnalytic(.e), .msg, fixed = TRUE)
  expect_false(.r$pd)
  expect_warning(foceiCovAnalytic(.e), .msg, fixed = TRUE)
  expect_identical(.e$covMethod, "r,s (full)")
  expect_equal(.e$cov, .anBad[1:2, 1:2])
})

test_that(".foceiFitInteraction() takes the interaction of the fit's likelihood", {
  # saem, nlme, vae and vi fits: the FOCEi control has the method's likelihood, while the
  # control their finalUi keeps is the output step's (interaction = 0)
  .ui <- new.env(parent = emptyenv())
  .ui$control <- list(interaction = 0L)
  expect_identical(.foceiFitInteraction(list(foceiControl = list(interaction = 1L)), .ui), 1L)
  expect_identical(.foceiFitInteraction(list(foceiControl = list(interaction = 0L)), .ui), 0L)
  # no usable FOCEi control: the ui's
  expect_identical(.foceiFitInteraction(list(), .ui), 0L)
  expect_identical(.foceiFitInteraction(list(foceiControl = list(interaction = NA_integer_)), .ui), 0L)
  expect_identical(.foceiFitInteraction(list(foceiControl = list(interaction = 1:2)), .ui), 0L)
  # a FOCEi control that cannot be built
  .err <- new.env(parent = emptyenv())
  makeActiveBinding("foceiControl", function() stop("no control"), .err)
  expect_identical(.foceiFitInteraction(.err, .ui), 0L)
  # and FOCEI when neither says
  .ui$control <- list()
  expect_identical(.foceiFitInteraction(list(), .ui), 1L)
})
