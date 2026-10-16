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
  expect_identical(.w, "\"linFim\" is not positive definite; installed as \"|linFim|\" (sqrtm)")
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
  expect_identical(.w, "\"linFim\" not usable; trying \"Ha\"")
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
  expect_identical(.w, "\"linFim\" not usable; trying \"Ha\"")
  expect_identical(.e$covMethod, "Ha")
  expect_equal(unname(.e$cov), unname(solve(.saemHa)))
})

test_that("a linFim covariance of the wrong size does not lend its variance block to the fallback", {
  .e <- .saemCovEnv()
  .lin <- diag(3) * 0.01
  attr(.lin, "varCov") <- matrix(0.01, 1, 1, dimnames = list("add.sd", "add.sd"))
  local_mocked_bindings(calc.COV = function(x) .lin)
  .w <- suppressMessages(capture_warnings(.saemCalcCov(.e)))
  expect_identical(.w, "\"linFim\" not usable; trying \"Ha\"")
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
      "\"linFim\" failed; trying \"Ha\"",
      "\"Ha\" is not finite; none installed"
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
  expect_identical(
    capture_warnings(.saemCalcCov(.e)),
    "\"Ha\" is not positive definite; installed as \"|Ha|\" (sqrtm)"
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
    "\"Ha\" has no estimated mu-referenced theta; none installed"
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
  expect_identical(
    capture_warnings(.saemCalcCov(.e)),
    "no \"Ha\" SE for non-mu-referenced theta(s) tka"
  )
  .ref <- solve(.saemHaKernel[1:2, 1:2])
  expect_identical(dimnames(.e$cov), list(c("tcl", "tv"), c("tcl", "tv")))
  expect_equal(.e$cov["tcl", "tcl"], .ref[1, 1])
  expect_equal(.e$cov["tcl", "tv"], .ref[1, 2])
  expect_equal(.e$cov["tv", "tv"], .ref[2, 2])
  expect_identical(.e$covMethod, "Ha")
})

test_that("the SE-less theta list is clipped to one line, the label is not", {
  .pre <- "no \"Ha\" SE for non-mu-referenced theta(s) "
  .r <- list(cov = .saemHa, label = "Ha", phi0 = strrep("t", 80))
  dimnames(.r$cov) <- list(.saemTheta, .saemTheta)
  .w <- .covSite(".saemInstallThetaCov", capture_warnings(.saemInstallThetaCov(.saemCovEnv(), .r)))
  expect_identical(.w, paste0(.pre, strrep("t", 29), "..."))
  expect_identical(cli::ansi_nchar(.w, type = "width"), 74L)
  .r$phi0 <- paste0("theta", 1:30)
  .w <- capture_warnings(.saemInstallThetaCov(.saemCovEnv(), .r))
  expect_identical(.w, paste0(.pre, "theta1, theta2, theta3, +27 more"))
})

test_that("each label a SAEM theta covariance carries gives whole warning lines", {
  .r <- list(cov = .saemHa, phi0 = "tka")
  dimnames(.r$cov) <- list(.saemTheta, .saemTheta)
  for (.l in c("Ha", "|Ha|", "linFim", "|linFim|")) {
    .r$label <- .l
    .w <- capture_warnings(.saemInstallThetaCov(.saemCovEnv(), .r))
    .want <- sprintf("no \"%s\" SE for non-mu-referenced theta(s) tka", .l)
    if (!identical(.covUnrepaired(.l), .l)) {
      .want <- c(.want, sprintf("\"%s\" is not positive definite; installed as \"%s\" (sqrtm)", .covUnrepaired(.l), .l))
    }
    expect_identical(.w, .want, label = .l)
    expect_true(.covSiteWhole(.w), label = .l)
  }
})

test_that("the linFim warning is clipped to one line when the next method is long", {
  .w <- .covSite(".saemLinFimWarn", capture_warnings(.saemLinFimWarn("not usable", strrep("m", 80))))
  .full <- paste0("\"linFim\" not usable; trying \"", strrep("m", 80), "\"")
  expect_identical(.w, paste0(substr(.full, 1L, 71L), "..."))
  expect_identical(cli::ansi_nchar(.w, type = "width"), 74L)
})

test_that("the linFim warning is one whole line for every method it can fall back to", {
  .next <- .covRealLabels(decorated = FALSE, full = FALSE)
  for (.why in c("failed", "not usable")) {
    .w <- vapply(.next, function(m) capture_warnings(.saemLinFimWarn(.why, m)), character(1), USE.NAMES = FALSE)
    expect_identical(.w, sprintf("\"linFim\" %s; trying \"%s\"", .why, .next))
    expect_true(.covSiteWhole(.w))
  }
  expect_identical(capture_warnings(.saemLinFimWarn("failed", NA_character_)), "\"linFim\" failed; none installed")
})

test_that("no population parameters is a covariance warning of one line", {
  .e <- .saemCovEnv(theta = character(0))
  .w <- .covSite(".saemCovRung", capture_warnings(.saemCalcCov(.e)))
  expect_identical(.w, "no population parameters in the model, no covariance matrix calculated")
  expect_lte(cli::ansi_nchar(.w, type = "width"), 74L)
  expect_identical(.e$covMethod, "none")
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
    "\"Ha\" has information rows not matching the thetas; none installed",
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
    "\"linFim (full)\" is not finite; kept \"linFim\"",
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

test_that("the SAEM covariance falls back only as covFallback lists", {
  # linFim fails; with no fallback listed there is no covariance
  .e <- .saemCovEnv()
  .e$ui$control$covFallback <- list(linFim = character(0))
  local_mocked_bindings(calc.COV = function(x) stop("ill-conditioned"))
  .w <- capture_warnings(.saemCalcCov(.e))
  expect_identical(.w, "\"linFim\" failed; none installed")
  expect_null(.e$cov)
  expect_null(.e$covMethod)
  expect_identical(.e$covTried, data.frame(method = "linFim", outcome = "not usable"))
  # with the default, the information matrix
  .e <- .saemCovEnv()
  .w <- capture_warnings(.saemCalcCov(.e))
  expect_identical(.w, "\"linFim\" failed; trying \"Ha\"")
  expect_identical(.e$covMethod, "Ha")
  expect_identical(.e$covTried, data.frame(method = c("linFim", "Ha"), outcome = c("not usable", "used")))
  # "r,s" (the information matrix) has no fallback by default
  expect_null(saemControl()$covFallback[["r,s"]])
})

test_that(".saemPas() gives one gain per iteration, a phase of 0 iterations none", {
  # burn-in at gain 1, then 1/k continuing past the burn-in's last k
  expect_equal(.saemPas(c(2L, 3L), 0:1), c(1, 1, 1 / 2, 1 / 3, 1 / 4))
  expect_equal(.saemPas(c(0L, 3L), 0:1), c(1, 1 / 2, 1 / 3))
  expect_equal(.saemPas(c(2L, 0L), 0:1), c(1, 1))
  expect_identical(.saemPas(c(0L, 0L), 0:1), numeric(0))
})
