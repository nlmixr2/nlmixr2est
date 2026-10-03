# The SAEM covariance step and its installers (R/saem.R), without a fit: the
# chol -> sqrtm repair, the fallbacks it reports, and the full-matrix install.

.saemTheta <- c("tka", "tcl")
.saemHa <- matrix(c(40, 5, 5, 30), 2) # SAEM information, theta block

.saemCovEnv <- function(
  Ha = .saemHa,
  covMethod = "linFim",
  theta = .saemTheta,
  fixed = rep(FALSE, length(theta)),
  cfg = NULL
) {
  .e <- new.env(parent = emptyenv())
  .e$ui <- new.env(parent = emptyenv())
  .e$ui$control <- list(covMethod = covMethod, covFull = TRUE)
  .e$ui$saemParamsToEstimate <- theta
  .e$ui$saemFixed <- fixed
  .e$ui$nonMuEtas <- character(0)
  .e$ui$iniDf <- data.frame(
    name = c(theta, "add.sd"),
    ntheta = seq_len(length(theta) + 1L),
    err = c(rep(NA, length(theta)), "add"),
    fix = c(fixed, FALSE)
  )
  .e$ui$mixProbs <- character(0)
  .e$saem <- list(Ha = Ha)
  attr(.e$saem, "saem.cfg") <- cfg
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
  expect_identical(.e$covMethod, "Ha")
  expect_false(exists(".saemFullCov", envir = .e, inherits = FALSE))
})

test_that("a linFim covariance that is not symmetric is not installed", {
  .e <- .saemCovEnv()
  .lin <- matrix(c(0.02, 0.001, 0.003, 0.01), 2) # chol() reads one triangle only
  attr(.lin, "varCov") <- matrix(0.01, 1, 1, dimnames = list("add.sd", "add.sd"))
  local_mocked_bindings(calc.COV = function(x) .lin)
  .w <- suppressMessages(capture_warnings(.saemCalcCov(.e)))
  expect_identical(.w, "linearization of FIM could not be used to calculate covariance")
  expect_identical(.e$covMethod, "Ha")
  expect_equal(unname(.e$cov), unname(solve(.saemHa)))
})

test_that("a linFim covariance of the wrong size does not lend its variance block to the fallback", {
  .e <- .saemCovEnv()
  .lin <- diag(3) * 0.01
  attr(.lin, "varCov") <- matrix(0.01, 1, 1, dimnames = list("add.sd", "add.sd"))
  local_mocked_bindings(calc.COV = function(x) .lin)
  .w <- suppressMessages(capture_warnings(.saemCalcCov(.e)))
  expect_identical(.w, "linearization of FIM could not be used to calculate covariance")
  expect_identical(.e$covMethod, "Ha")
  expect_equal(unname(.e$cov), unname(solve(.saemHa)))
  expect_false(exists(".saemFullCov", envir = .e, inherits = FALSE))
  expect_false(exists(".saemCovMethod", envir = .e, inherits = FALSE))
})

test_that("an unusable information matrix is reported instead of silently leaving no covariance", {
  .e <- .saemCovEnv(Ha = matrix(c(1, NaN, NaN, 1), 2))
  local_mocked_bindings(calc.COV = function(x) stop("singular"))
  .w <- suppressMessages(capture_warnings(.saemCalcCov(.e)))
  expect_identical(
    .w,
    c(
      "SAEM covariance by linearization failed; using the SAEM information matrix",
      "\"Ha\" covariance is not finite; none installed"
    )
  )
  expect_false(exists("cov", envir = .e, inherits = FALSE))
  expect_false(exists("covMethod", envir = .e, inherits = FALSE))
})

test_that("covMethod r,s/r/s install the inverse of Ha's theta block under its own name", {
  for (.cm in c("r,s", "r", "s")) {
    .e <- .saemCovEnv(covMethod = .cm)
    expect_silent(.saemCalcCov(.e))
    expect_equal(.e$cov, solve(.saemHa), ignore_attr = TRUE)
    expect_identical(dimnames(.e$cov), list(.saemTheta, .saemTheta))
    expect_identical(.e$covMethod, "Ha")
  }
})

test_that("an indefinite Ha theta block is repaired by sqrtm and labelled so", {
  .bad <- matrix(c(1, 2, 2, 1), 2)
  .e <- .saemCovEnv(Ha = .bad, covMethod = "r,s")
  expect_warning(
    .saemCalcCov(.e),
    "covariance matrix non-positive definite, corrected by sqrtm(Ha %*% Ha)",
    fixed = TRUE
  )
  expect_identical(.e$covMethod, "|Ha|")
  expect_equal(unname(.e$cov), unname(solve(sqrtm(.bad %*% .bad))))
})

# Ha laid out the way src/saem.cpp lays it out: the structural rows first in
# [phi1 mu][phi0 mu] order (a row for a fixed theta too), then log-Omega and
# log-sigma2 rows.  .saemHaKernel is the matching information in that order.
.saemHaKernel <- rbind(
  c(40, 5, 2, 1, 0.5),
  c(5, 30, 3, 0.4, 0.2),
  c(2, 3, 20, 0.3, 0.1),
  c(1, 0.4, 0.3, 10, 0.6),
  c(0.5, 0.2, 0.1, 0.6, 8)
)

test_that("the Ha theta block skips a fixed theta's row and names rows in kernel order", {
  # kernel order [phi1 mu] = tka, tcl, tv (no phi0); tcl is fixed
  .e <- .saemCovEnv(
    Ha = .saemHaKernel,
    covMethod = "r,s",
    theta = c("tka", "tcl", "tv"),
    fixed = c(FALSE, TRUE, FALSE),
    cfg = list(i1 = 0:2, i0 = integer(0), nphi1 = 3L, nphi0 = 0L)
  )
  expect_silent(.saemCalcCov(.e))
  .ref <- solve(.saemHaKernel[c(1, 3), c(1, 3)])
  expect_identical(dimnames(.e$cov), list(c("tka", "tv"), c("tka", "tv")))
  expect_equal(.e$cov["tka", "tka"], .ref[1, 1])
  expect_equal(.e$cov["tka", "tv"], .ref[1, 2])
  expect_equal(.e$cov["tv", "tv"], .ref[2, 2])
})

test_that("a Ha theta block with one estimated theta is its 1 x 1 inverse", {
  .e <- .saemCovEnv(
    Ha = .saemHaKernel,
    covMethod = "r,s",
    theta = c("tka", "tcl", "tv"),
    fixed = c(TRUE, FALSE, TRUE),
    cfg = list(i1 = 0:2, i0 = integer(0), nphi1 = 3L, nphi0 = 0L)
  )
  expect_silent(.saemCalcCov(.e))
  expect_identical(.e$cov, matrix(1 / 30, 1, 1, dimnames = list("tcl", "tcl")))
  expect_identical(.e$covMethod, "Ha")
})

test_that("a Ha theta block with every theta fixed installs nothing and says so", {
  .e <- .saemCovEnv(
    Ha = .saemHaKernel,
    covMethod = "r,s",
    theta = c("tka", "tcl", "tv"),
    fixed = c(TRUE, TRUE, TRUE),
    cfg = list(i1 = 0:2, i0 = integer(0), nphi1 = 3L, nphi0 = 0L)
  )
  expect_identical(
    capture_warnings(.saemCalcCov(.e)),
    "\"Ha\" covariance could not be computed (no mu-referenced theta is estimated); none installed"
  )
  expect_false(exists("cov", envir = .e, inherits = FALSE))
  expect_false(exists("covMethod", envir = .e, inherits = FALSE))
})

test_that("the Ha theta block drops a non-mu-referenced theta and reports it", {
  # tka has no eta (phi0), so the kernel order is [phi1 mu][phi0 mu] = tcl, tv, tka
  .e <- .saemCovEnv(
    Ha = .saemHaKernel,
    covMethod = "r",
    theta = c("tka", "tcl", "tv"),
    cfg = list(i1 = 1:2, i0 = 0L, nphi1 = 2L, nphi0 = 1L)
  )
  expect_warning(
    .saemCalcCov(.e),
    "\"Ha\" covariance has no row for the non-mu-referenced theta(s) tka; they have no standard error",
    fixed = TRUE
  )
  .ref <- solve(.saemHaKernel[1:2, 1:2])
  expect_identical(dimnames(.e$cov), list(c("tcl", "tv"), c("tcl", "tv")))
  expect_equal(.e$cov["tcl", "tcl"], .ref[1, 1])
  expect_equal(.e$cov["tcl", "tv"], .ref[1, 2])
  expect_equal(.e$cov["tv", "tv"], .ref[2, 2])
  expect_identical(.e$covMethod, "Ha")
})

test_that("a Ha theta block whose rows cannot be ordered is not installed", {
  # a phi0 theta with no verified [phi1][phi0] partition (here: a covariate
  # coefficient lengthens the parameter list past the phi columns)
  .e <- .saemCovEnv(
    Ha = .saemHaKernel,
    covMethod = "s",
    theta = c("tka", "tcl", "wt.cl", "tv"),
    cfg = list(i1 = 1:2, i0 = 0L, nphi1 = 2L, nphi0 = 1L)
  )
  expect_warning(
    .saemCalcCov(.e),
    paste0(
      "\"Ha\" covariance could not be computed (the information rows cannot be ",
      "matched to the thetas); none installed"
    ),
    fixed = TRUE
  )
  expect_false(exists("cov", envir = .e, inherits = FALSE))
})

test_that("saemControl() reads an integer covMethod as a foceiControl() slot", {
  expect_identical(saemControl(covMethod = 0L)$covMethod, "")
  expect_identical(saemControl(covMethod = 1L)$covMethod, "r,s")
  expect_identical(saemControl(covMethod = 2L)$covMethod, "r")
  expect_identical(saemControl(covMethod = 3L)$covMethod, "s")
  expect_error(saemControl(covMethod = 5L), "foceiControl() slot", fixed = TRUE)
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
  # the sqrtm repair was reported when the matrix was computed
  expect_no_warning(.saemInstallFullCov(.e))
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
