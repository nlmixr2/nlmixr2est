# .foceiInstallFdFullCov routing unit tests (always-run core: no fit, just the R install
# assembling fit$cov / covR / covS / covRS from the stashed C++ pieces .fdFullCov (Rinv_full)
# and .fdFullS (Sfull) per the fit's covMethod string).

.mkFdEnv <- function(covMethod, Rinv, S = NULL) {
  .e <- new.env(parent = emptyenv())
  .nm <- c("tka", "om.eta.ka")
  dimnames(Rinv) <- list(.nm, .nm)
  .e$.fdFullCov <- Rinv
  if (!is.null(S)) { dimnames(S) <- list(.nm, .nm); .e$.fdFullS <- S }
  .e$covMethod <- covMethod
  .e
}

# well-conditioned PD pieces
.Rinv <- matrix(c(4, 1, 1, 3), 2)
.S    <- matrix(c(2, 0.5, 0.5, 1), 2)
.sandwich <- .Rinv %*% .S %*% .Rinv

test_that("covMethod='r,s' installs the true full sandwich Rinv %*% S %*% Rinv", {
  .e <- .mkFdEnv("r,s", .Rinv, .S)
  .foceiInstallFdFullCov(.e)
  expect_equal(unname(.e$cov), .sandwich)
  expect_equal(unname(.e$covR), unname(.Rinv))
  expect_equal(unname(.e$covS), unname(solve(.S)))
  expect_equal(unname(.e$covRS), .sandwich)
  expect_identical(rownames(.e$cov), c("tka", "om.eta.ka"))
})

test_that("sandwich variant covMethod strings route to r,s", {
  for (.cm in c("r+,s", "|r|,|s|", "r,s+", "r+,s+", "|r|,s")) {
    .e <- .mkFdEnv(.cm, .Rinv, .S)
    .foceiInstallFdFullCov(.e)
    expect_equal(unname(.e$cov), .sandwich, info = .cm)
  }
})

test_that("covMethod='s' installs solve(Sfull)", {
  .e <- .mkFdEnv("s", .Rinv, .S)
  .foceiInstallFdFullCov(.e)
  expect_equal(unname(.e$cov), unname(solve(.S)))
  expect_equal(unname(.e$covS), unname(solve(.S)))
})

test_that("covMethod='r' installs Rinv and does not touch covS/covRS", {
  .e <- .mkFdEnv("r", .Rinv)          # no S stashed
  .foceiInstallFdFullCov(.e)
  expect_equal(unname(.e$cov), unname(.Rinv))
  expect_equal(unname(.e$covR), unname(.Rinv))
  expect_false(exists("covS", envir = .e, inherits = FALSE))
  expect_false(exists("covRS", envir = .e, inherits = FALSE))
})

test_that("non-FD covMethod (analytic/failed/boundary/empty) is a no-op", {
  for (.cm in c("analytic", "failed", "", "Boundary issue")) {
    .e <- .mkFdEnv(.cm, .Rinv, .S)
    .foceiInstallFdFullCov(.e)
    expect_false(exists("cov", envir = .e, inherits = FALSE), info = .cm)
  }
})

test_that("s/r,s with a missing or non-finite Sfull is a no-op", {
  .e <- .mkFdEnv("r,s", .Rinv)        # no .fdFullS
  .foceiInstallFdFullCov(.e)
  expect_false(exists("cov", envir = .e, inherits = FALSE))
  .Sbad <- matrix(c(1, NA, NA, 1), 2)
  .e2 <- .mkFdEnv("s", .Rinv, .Sbad)
  .foceiInstallFdFullCov(.e2)
  expect_false(exists("cov", envir = .e2, inherits = FALSE))
})

test_that("PD guard rejects an indefinite assembled cov", {
  .Rbad <- matrix(c(1, 0, 0, -2), 2)  # negative variance -> not PD
  .e <- .mkFdEnv("r", .Rbad)
  .foceiInstallFdFullCov(.e)
  expect_false(exists("cov", envir = .e, inherits = FALSE))
})

test_that("no .fdFullCov stashed is a no-op", {
  .e <- new.env(parent = emptyenv())
  .e$covMethod <- "r,s"
  .foceiInstallFdFullCov(.e)
  expect_false(exists("cov", envir = .e, inherits = FALSE))
})

# --- covariance scope naming + the covList cache the install populates --------

test_that("covariance scope names round-trip", {
  expect_identical(.covFullName("s"), "s (full)")
  expect_identical(.covFullName("r,s"), "r,s (full)")
  expect_identical(.covFullName("analytic"), "analytic (full)")
  expect_identical(.covFullName("s (full)"), "s (full)")       # idempotent
  expect_identical(.covBaseName("r,s (full)"), "r,s")
  expect_identical(.covBaseName("r,s"), "r,s")
  expect_true(.covIsFull("analytic (full)"))
  expect_false(.covIsFull("analytic"))
  expect_false(.covIsFull(character(0)))
  expect_false(.covIsFull(NA_character_))
})

test_that("decorated native covMethod strings resolve to their FD type", {
  expect_identical(.covFdType("|r|,s+"), "r,s")
  expect_identical(.covFdType("|s| (full)"), "s")
  expect_identical(.covFdType("r+"), "r")
  expect_identical(.covFdType("analytic"), "")
  expect_identical(.covFdType("failed"), "")
})

test_that(".covSameName ignores decorations but not the scope", {
  expect_true(.covSameName("s", "|s|"))
  expect_true(.covSameName("r,s (full)", "r+,s (full)"))
  expect_false(.covSameName("s", "s (full)"))
  expect_false(.covSameName("r", "s"))
})

test_that("the FD-full install names the full shape and caches the others", {
  .e <- .mkFdEnv("r,s", .Rinv, .S)
  expect_true(.foceiInstallFdFullCov(.e))
  expect_identical(.e$covMethod, "r,s (full)")
  # the installed shape is never left in the cache; the other full shapes are
  expect_false("r,s (full)" %in% names(.e$covList))
  expect_setequal(names(.e$covList), c("r (full)", "s (full)"))
  expect_equal(unname(.e$covList[["s (full)"]]), unname(solve(.S)))
})

test_that("the native theta-only covariance is cached under its own type", {
  # covMethod="s" writes only e$cov in C++ (no e$covS), so the installed native
  # covariance is the only copy of the theta-only shape -- cache it before the swap
  .e <- .mkFdEnv("s", .Rinv, .S)
  .native <- matrix(c(0.25, 0.01, 0.01, 0.16), 2, dimnames = list(c("tka", "om.eta.ka"),
                                                                  c("tka", "om.eta.ka")))
  .e$cov <- .native
  expect_true(.foceiInstallFdFullCov(.e))
  expect_identical(.e$covMethod, "s (full)")
  expect_equal(.e$covList[["s"]], .native)
  expect_false("s (full)" %in% names(.e$covList))
})

test_that("a decorated native covMethod keeps its decoration under the suffix", {
  .e <- .mkFdEnv("|s|", .Rinv, .S)
  expect_true(.foceiInstallFdFullCov(.e))
  expect_identical(.e$covMethod, "|s| (full)")
})

test_that("a no-op install caches nothing", {
  .e <- .mkFdEnv("analytic", .Rinv, .S)
  expect_false(.foceiInstallFdFullCov(.e))
  expect_false(exists("covList", envir = .e, inherits = FALSE))
})
