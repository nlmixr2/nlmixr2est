# .foceiInstallFdFullCov() checks every full-shape covariance it stores or caches
# (R/foceiCovFdFull.R); no fit, just the R install from stashed C++ pieces.

.fdNames <- c("tka", "om.eta.ka")

.fdEnv <- function(covMethod, Rinv, S = NULL) {
  .e <- new.env(parent = emptyenv())
  dimnames(Rinv) <- list(.fdNames, .fdNames)
  .e$.fdFullCov <- Rinv
  if (!is.null(S)) {
    dimnames(S) <- list(.fdNames, .fdNames)
    .e$.fdFullS <- S
  }
  # what the native theta-only step installed
  .e$cov <- matrix(0.25, 1, 1, dimnames = list("tka", "tka"))
  .e$covR <- matrix(0.2, 1, 1, dimnames = list("tka", "tka"))
  .e$covMethod <- covMethod
  .e
}

.fdS <- matrix(c(2, 0.5, 0.5, 1), 2)
.fdRbad <- matrix(c(4, 1, 1, -0.5), 2) # indefinite

test_that(".foceiFdFullShapes() judges the sandwich by its pieces", {
  .sh <- .foceiFdFullShapes(.fdRbad, .fdS)
  expect_false(.sh$r$ok)
  expect_identical(.sh$r$reason, "is not positive definite")
  expect_true(.sh$s$ok)
  # R^-1 S R^-1 is positive definite for any non-singular R, so on its own it passes
  expect_true(.covGuard(.fdRbad %*% .fdS %*% .fdRbad)$ok)
  expect_false(.sh[["r,s"]]$ok)
  expect_identical(.sh[["r,s"]]$reason, "needs a positive-definite R")
  .sh <- .foceiFdFullShapes(matrix(c(4, 1, 1, 3), 2), matrix(c(1, 1, 1, 1), 2))
  expect_true(.sh$r$ok)
  expect_identical(.sh$s$reason, "could not be computed")
  expect_identical(.sh[["r,s"]]$reason, "needs a positive-definite S")
  .sh <- .foceiFdFullShapes(matrix(c(4, 1, 1, 3), 2), NULL)
  expect_true(.sh$r$ok)
  expect_false(.sh$s$ok)
})

test_that("an indefinite full R is not installed inside a positive-definite sandwich", {
  .e <- .fdEnv("r,s", .fdRbad, .fdS)
  expect_warning(
    .ok <- .foceiInstallFdFullCov(.e),
    "\"r,s (full)\" covariance needs a positive-definite R; kept \"r,s\"",
    fixed = TRUE
  )
  expect_false(.ok)
  expect_identical(.e$covMethod, "r,s")
  expect_equal(.e$cov, matrix(0.25, 1, 1, dimnames = list("tka", "tka")))
  expect_equal(.e$covR, matrix(0.2, 1, 1, dimnames = list("tka", "tka")))
  # the usable shapes stay swappable: the native R and the full S
  expect_identical(sort(names(.e$covList)), c("r", "s (full)"))
  expect_equal(unname(.e$covList[["s (full)"]]), unname(solve(.fdS)))
})

test_that("an indefinite full R stays out of covR and covList when another shape installs", {
  .e <- .fdEnv("s", .fdRbad, .fdS)
  expect_true(.foceiInstallFdFullCov(.e))
  expect_identical(.e$covMethod, "s (full)")
  expect_equal(unname(.e$cov), unname(solve(.fdS)))
  expect_equal(unname(.e$covS), unname(solve(.fdS)))
  expect_equal(.e$covR, matrix(0.2, 1, 1, dimnames = list("tka", "tka")))
  expect_false(exists("covRS", envir = .e, inherits = FALSE))
  # the native shapes are cached; of the full ones only the usable S is (and it is installed)
  expect_identical(sort(names(.e$covList)), c("r", "s"))
})

test_that("the installed full covariance refreshes the condition numbers", {
  .R <- matrix(c(4, 1, 1, 3), 2)
  .e <- .fdEnv("r", .R)
  .e$objDf <- data.frame(OBJF = 1, "Condition#(Cov)" = 99, "Condition#(Cor)" = 99, check.names = FALSE)
  expect_true(.foceiInstallFdFullCov(.e))
  .ev <- eigen(.R, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(.e$objDf[["Condition#(Cov)"]], max(.ev) / min(.ev))
  expect_equal(.e$conditionNumberCov, max(.ev) / min(.ev))
  expect_equal(unname(.e$fullCor), unname(.cov2cor(.R)))
})
