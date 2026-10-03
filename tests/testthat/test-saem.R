# The SAEM covariance step and its installers (R/saem.R), without a fit: the
# chol -> sqrtm repair, the fallbacks it reports, and the full-matrix install.

.saemTheta <- c("tka", "tcl")
.saemHa <- matrix(c(40, 5, 5, 30), 2) # SAEM information, theta block

.saemCovEnv <- function(Ha = .saemHa, covMethod = "linFim") {
  .e <- new.env(parent = emptyenv())
  .e$ui <- new.env(parent = emptyenv())
  .e$ui$control <- list(covMethod = covMethod, covFull = TRUE)
  .e$ui$saemParamsToEstimate <- .saemTheta
  .e$ui$saemFixed <- c(FALSE, FALSE)
  .e$ui$iniDf <- data.frame(
    name = c("tka", "tcl", "add.sd"),
    ntheta = c(1L, 2L, 3L),
    err = c(NA, NA, "add"),
    fix = c(FALSE, FALSE, FALSE)
  )
  .e$ui$mixProbs <- character(0)
  .e$saem <- list(Ha = Ha)
  .e
}

test_that(".saemCovRepair() factors, repairs or gives up", {
  .r <- .saemCovRepair(.saemHa)
  expect_false(.r$sqrtm)
  expect_identical(.r$mat, .saemHa)
  .bad <- matrix(c(1, 2, 2, 1), 2)
  .r <- .saemCovRepair(.bad)
  expect_true(.r$sqrtm)
  expect_equal(.r$mat, sqrtm(.bad %*% t(.bad)))
  expect_equal(eigen(.r$mat, symmetric = TRUE, only.values = TRUE)$values, c(3, 1))
  expect_null(suppressMessages(.saemCovRepair(matrix(c(1, NaN, NaN, 1), 2))))
  # a partial factorization judges only the identified block
  .na <- matrix(c(1, NA, NA, NA), 2)
  expect_identical(.saemCovRepair(.na, partial = TRUE)$mat, .na)
})

test_that("a linFim covariance repaired by sqrtm keeps that label with its full matrix", {
  .e <- .saemCovEnv()
  .lin <- matrix(c(0.02, 0.03, 0.03, 0.01), 2) # indefinite
  attr(.lin, "varCov") <- matrix(0.01, 1, 1, dimnames = list("add.sd", "add.sd"))
  local_mocked_bindings(calc.COV = function(x) .lin)
  .w <- capture_warnings(.saemCalcCov(.e))
  expect_identical(.w, "covariance matrix non-positive definite, corrected by sqrtm(linFim %*% linFim)")
  expect_identical(.e$covMethod, "|linFim|")
  .abs <- sqrtm(.lin %*% t(.lin))
  expect_equal(unname(.e$cov), unname(.abs))
  expect_identical(.e$.saemCovMethod, "|linFim|")
  expect_equal(unname(.e$.saemFullCov), unname(rbind(cbind(.abs, 0), c(0, 0, 0.01))))
})

test_that("a linFim covariance that cannot be used is reported, and the information matrix used", {
  .e <- .saemCovEnv()
  local_mocked_bindings(calc.COV = function(x) matrix(c(1, NaN, NaN, 1), 2))
  .w <- suppressMessages(capture_warnings(.saemCalcCov(.e)))
  expect_identical(.w, "linearization of FIM could not be used to calculate covariance")
  expect_equal(unname(.e$cov), unname(solve(.saemHa)))
  expect_identical(dimnames(.e$cov), list(.saemTheta, .saemTheta))
  expect_false(exists("covMethod", envir = .e, inherits = FALSE))
  expect_false(exists(".saemFullCov", envir = .e, inherits = FALSE))
})

test_that("an unusable information matrix is reported instead of silently leaving no covariance", {
  .e <- .saemCovEnv(Ha = matrix(c(1, NaN, NaN, 1), 2))
  local_mocked_bindings(calc.COV = function(x) stop("singular"))
  .w <- suppressMessages(capture_warnings(.saemCalcCov(.e)))
  expect_identical(
    .w,
    c(
      "SAEM covariance by linearization failed; using the SAEM information matrix",
      "FIM non-positive definite and cannot be used to calculate the covariance"
    )
  )
  expect_false(exists("cov", envir = .e, inherits = FALSE))
  expect_false(exists("covMethod", envir = .e, inherits = FALSE))
})

test_that("a non-linFim covMethod inverts the information matrix and sets no label", {
  .e <- .saemCovEnv(covMethod = "r,s")
  expect_silent(.saemCalcCov(.e))
  expect_equal(unname(.e$cov), unname(solve(.saemHa)))
  expect_false(exists("covMethod", envir = .e, inherits = FALSE))
})

.saemFullEnv <- function(full, label = "|linFim|") {
  .e <- new.env(parent = emptyenv())
  .e$.saemFullCov <- full
  .e$.saemCovMethod <- label
  .e$cov <- full[1:2, 1:2]
  .e$covMethod <- "failed" # a stale label left by the shared finalization
  .e$parFixedDf <- data.frame(
    Estimate = c(0.45, 1, 0.7),
    SE = c(0.5, 0.5, NA_real_),
    "%RSE" = c(NA_real_, NA_real_, NA_real_),
    check.names = FALSE,
    row.names = c("tka", "tcl", "add.sd")
  )
  .e$objDf <- data.frame(OBJF = 1, "Condition#(Cov)" = 99, "Condition#(Cor)" = 99, check.names = FALSE)
  .e
}

.saemFull <- matrix(
  c(0.04, 0.01, 0, 0.01, 0.09, 0, 0, 0, 0.0025),
  3,
  dimnames = list(c("tka", "tcl", "add.sd"), c("tka", "tcl", "add.sd"))
)

test_that("the full SAEM covariance installs under the label it was computed with", {
  .e <- .saemFullEnv(.saemFull)
  .saemInstallFullCov(.e)
  expect_identical(.e$covMethod, "|linFim|")
  expect_equal(.e$cov, .saemFull)
  # only the SE the theta-only table lacked is filled
  expect_equal(.e$parFixedDf$SE, c(0.5, 0.5, 0.05))
  .ev <- eigen(.saemFull, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(.e$objDf[["Condition#(Cov)"]], max(.ev) / min(.ev))
  expect_equal(.e$conditionNumberCov, max(.ev) / min(.ev))
})

test_that("a full SAEM covariance that is not usable is reported and the theta block kept", {
  .bad <- .saemFull
  .bad[3, 3] <- NA
  .e <- .saemFullEnv(.bad, "linFim")
  .e$covMethod <- "linFim"
  expect_warning(
    .saemInstallFullCov(.e),
    "\"linFim (full)\" covariance is not finite; kept \"linFim\"",
    fixed = TRUE
  )
  expect_equal(.e$cov, .saemFull[1:2, 1:2])
  expect_identical(.e$objDf[["Condition#(Cov)"]], 99)
})

test_that("the saem finalization control substitutes fixed parameters as the fit did", {
  .finalControl <- function(ctl) {
    .e <- new.env(parent = emptyenv())
    .e$saemControl <- ctl
    .e$.etaMat <- matrix(0, 2, 1)
    .e$ui <- list(foceiSkipCov = NULL)
    .saemControlToFoceiControl(.e, assign = FALSE)
  }
  for (.lf in c(TRUE, FALSE)) {
    .fc <- .finalControl(saemControl(literalFix = .lf))
    expect_identical(.fc$literalFix, .lf)
    # saem never substitutes the fixed residual parameters
    expect_false(.fc$literalFixRes)
  }
  # the default: saem keeps fixed thetas in the model
  expect_false(.finalControl(saemControl())$literalFix)
  # a control saved before saemControl() had literalFix
  .old <- saemControl()
  .old$literalFix <- NULL
  expect_false(.finalControl(.old)$literalFix)
})
