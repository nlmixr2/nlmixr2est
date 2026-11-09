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

test_that("an indefinite full R falls back from the sandwich to the full S (#1152)", {
  .e <- .fdEnv("r,s", .fdRbad, .fdS)
  .w <- .covSite(".foceiInstallFdFullCov", capture_warnings(.ok <- .foceiInstallFdFullCov(.e)))
  expect_identical(.w, "\"r,s (full)\" needs a positive-definite R; using \"s (full)\"")
  expect_true(.ok)
  expect_identical(.e$covMethod, "s (full)")
  expect_equal(unname(.e$cov), unname(solve(.fdS)))
  expect_equal(unname(.e$covS), unname(solve(.fdS)))
  expect_false(exists("covR", envir = .e, inherits = FALSE))
  expect_false(exists("covRS", envir = .e, inherits = FALSE))
  # the native theta-only shapes stay swappable; the unusable full ones are not cached
  expect_identical(sort(names(.e$covList)), c("r", "r,s"))
  expect_equal(.e$covList[["r"]], matrix(0.2, 1, 1, dimnames = list("tka", "tka")))
  expect_equal(.e$covList[["r,s"]], matrix(0.25, 1, 1, dimnames = list("tka", "tka")))
})

test_that("an indefinite full R stays out of covR and covList when another shape installs", {
  .e <- .fdEnv("s", .fdRbad, .fdS)
  expect_true(.foceiInstallFdFullCov(.e))
  expect_identical(.e$covMethod, "s (full)")
  expect_equal(unname(.e$cov), unname(solve(.fdS)))
  expect_equal(unname(.e$covS), unname(solve(.fdS)))
  expect_false(exists("covR", envir = .e, inherits = FALSE))
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

test_that("same() copies move with their master in the full FD stage (#1316)", {
  # eta.cl, eta.v, then three occasions of a correlated iov.cl/iov.v block
  .onm <- c("eta.cl", "eta.v", paste0("rx.", c("iov.cl", "iov.v"), ".", rep(1:3, each = 2)))
  .same <- c(0L, 0L, 0L, 0L, 3L, 4L, 3L, 4L)
  .pairs <- .omegaDeclaredPairs(
    diag(8),
    data.frame(
      name = c(.onm[1:2], .onm[3], "(c1)", .onm[4], .onm[5], "(c2)", .onm[6], .onm[7], "(c3)", .onm[8]),
      neta1 = c(1, 2, 3, 4, 4, 5, 6, 6, 7, 8, 8),
      neta2 = c(1, 2, 3, 3, 4, 5, 5, 6, 7, 7, 8),
      fix = FALSE
    )
  )
  .tie <- .foceiFdFullSameTie(.pairs, .same, .onm)
  expect_identical(sum(.tie$free), 5L)
  expect_identical(.tie$cpA, c(5L, 6L, 6L, 7L, 8L, 8L))
  expect_identical(.tie$cpB, c(5L, 5L, 6L, 7L, 7L, 8L))
  # each copy points at its master element among the free pairs
  .free <- .pairs[.tie$free, , drop = FALSE]
  expect_identical(.free[.tie$cpOf + 1L, , drop = FALSE], .pairs[c(3, 4, 5, 3, 4, 5), , drop = FALSE])
  expect_identical(
    .foceiOmegaCovNames(.free, .tie$onm),
    c("om.eta.cl", "om.eta.v", "om.iov.cl", "cov.iov.v.iov.cl", "om.iov.v")
  )
  # no same() map: untouched
  expect_true(all(.foceiFdFullSameTie(.pairs, NULL, .onm)$free))
})

test_that("a same() group is named by the parameter it repeats", {
  .old <- .uiIovEnv$iovMaster
  on.exit(.uiIovEnv$iovMaster <- .old)
  .uiIovEnv$iovMaster <- list()
  # occasion levels sharing a prefix
  expect_identical(.foceiSameGroupName("rx.iov.cl.1", c("rx.iov.cl.10", "rx.iov.cl.2")), "iov.cl")
  # a user's own same() block keeps its master's name
  expect_identical(.foceiSameGroupName("eta.cl", "eta.v"), "eta.cl")
  # the expansion's record wins over the name
  .uiIovEnv$iovMaster <- list(iov.cl = "rx.iov.cl.a.b")
  expect_identical(.foceiSameGroupName("rx.iov.cl.a.b", "rx.iov.cl.a.c"), "iov.cl")
})

test_that("a chained same() map and a fixed master tie as expected (#1316)", {
  .pairs <- cbind(1:3, 1:3)
  .onm <- c("eta.a", "eta.b", "eta.c")
  # eta.c repeats eta.b, which repeats eta.a
  .tie <- .foceiFdFullSameTie(.pairs, c(0L, 1L, 2L), .onm)
  expect_identical(.tie$free, c(TRUE, FALSE, FALSE))
  expect_identical(.tie$cpA, 2:3)
  expect_identical(.tie$cpOf, c(0L, 0L))
  # the master is fixed (not among the pairs): its copy is not a parameter and does not move
  .tie <- .foceiFdFullSameTie(cbind(c(2L, 3L), c(2L, 3L)), c(0L, 0L, 1L), .onm)
  expect_identical(.tie$free, c(TRUE, FALSE))
  expect_length(.tie$cpA, 0L)
})

test_that(".foceiFdFullParams() ties a user's same() blocks (#1316)", {
  .f <- function() {
    ini({
      tcl <- 1
      add.sd <- 1
      eta.a + eta.b ~ c(0.3, 0.1, 0.2)
      eta.c + eta.d ~ same()
      eta.e + eta.f ~ same()
    })
    model({
      cl <- tcl * exp(eta.a + eta.b + eta.c + eta.d + eta.e + eta.f)
      cl ~ add(add.sd)
    })
  }
  .ui <- suppressWarnings(suppressMessages(rxode2::rxode2(.f)))
  .e <- new.env(parent = emptyenv())
  .e$ui <- .ui
  .e$omega <- .ui$omega
  .p <- .foceiFdFullParams(.e)
  expect_identical(.p$names, c("tcl", "add.sd", "om.eta.a", "cov.eta.b.eta.a", "om.eta.b"))
  expect_identical(.p$omA, c(1L, 2L, 2L))
  expect_identical(.p$omB, c(1L, 1L, 2L))
  expect_identical(.p$cpA, c(3L, 4L, 4L, 5L, 6L, 6L))
  expect_identical(.p$cpB, c(3L, 3L, 4L, 5L, 5L, 6L))
  expect_identical(.p$cpOf, c(0L, 1L, 2L, 0L, 1L, 2L))
})
