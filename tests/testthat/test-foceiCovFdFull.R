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
  .tie <- .foceiFdFullSameTie(.pairs, rep(FALSE, nrow(.pairs)), .same, .onm, c("tcl", "iov.cl", "iov.v"))
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
  expect_true(all(.foceiFdFullSameTie(.pairs, rep(FALSE, nrow(.pairs)), NULL, .onm)$free))
})

test_that("a same() group is named by the parameter it repeats", {
  .th <- c("tcl", "iov.cl", "iov")
  # occasion levels sharing a prefix, or holding a dot
  expect_identical(.foceiSameGroupName("rx.iov.cl.1", c("rx.iov.cl.10", "rx.iov.cl.2"), .th), "iov.cl")
  expect_identical(.foceiSameGroupName("rx.iov.cl.day.1", "rx.iov.cl.day.2", .th), "iov.cl")
  # a user's own same() block keeps its master's name
  expect_identical(.foceiSameGroupName("eta.cl", "eta.v", .th), "eta.cl")
  expect_identical(.foceiSameGroupName("rx.eta.cl", "rx.eta.v", .th), "rx.eta.cl")
  # occasion levels that read as another theta: the expansion's record decides
  .old <- .uiIovEnv$iovMaster
  on.exit(.uiIovEnv$iovMaster <- .old)
  .uiIovEnv$iovMaster <- list(iov = "rx.iov.cl.1")
  expect_identical(.foceiSameGroupName("rx.iov.cl.1", "rx.iov.cl.2", .th), "iov")
})

test_that("chained, fixed and malformed same() maps tie as expected (#1316)", {
  .pairs <- cbind(1:3, 1:3)
  .onm <- c("eta.a", "eta.b", "eta.c")
  .ok <- rep(FALSE, 3)
  # eta.c repeats eta.b, which repeats eta.a
  .tie <- .foceiFdFullSameTie(.pairs, .ok, c(0L, 1L, 2L), .onm)
  expect_identical(.tie$free, c(TRUE, FALSE, FALSE))
  expect_identical(.tie$cpA, 2:3)
  expect_identical(.tie$cpOf, c(0L, 0L))
  # a fixed master's copies stay put; a fixed copy of a free master moves with it
  .tie <- .foceiFdFullSameTie(.pairs, c(TRUE, FALSE, FALSE), c(0L, 0L, 1L), .onm)
  expect_identical(.tie$free, c(FALSE, TRUE, FALSE))
  expect_length(.tie$cpA, 0L)
  .tie <- .foceiFdFullSameTie(.pairs, c(FALSE, FALSE, TRUE), c(0L, 0L, 1L), .onm)
  expect_identical(.tie$free, c(TRUE, TRUE, FALSE))
  expect_identical(.tie$cpA, 3L)
  expect_identical(.tie$cpOf, 0L)
  # a map rxode2 would not produce leaves the matrix untied
  for (.s in list(c(0L, 5L, 0L), c(2L, 1L, 0L), c(1L, 0L, 0L), c(0L, -1L, 1L))) {
    expect_identical(.foceiFdFullSameTie(.pairs, .ok, .s, .onm)$free, !.ok)
  }
  # two copies of one eta correlated inside a block; a copy correlated with a non-copy
  .p2 <- rbind(c(1L, 1L), c(2L, 2L), c(3L, 2L), c(3L, 3L))
  expect_identical(.foceiFdFullSameTie(.p2, rep(FALSE, 4), c(0L, 1L, 1L), .onm)$free, rep(TRUE, 4))
  expect_identical(.foceiFdFullSameTie(.p2, rep(FALSE, 4), c(0L, 0L, 1L), .onm)$free, rep(TRUE, 4))
  # a copy covariance its (uncorrelated) master does not declare
  .p3 <- rbind(c(1L, 1L), c(2L, 2L), c(3L, 3L), c(4L, 3L), c(4L, 4L))
  expect_identical(.foceiFdFullSameTie(.p3, rep(FALSE, 5), c(0L, 0L, 1L, 2L), paste0("e", 1:4))$free, rep(TRUE, 5))
  # an element spanning two copy blocks; a copy block missing its master's covariance
  .p5 <- .omegaDeclaredPairs(
    diag(6),
    data.frame(
      name = paste0("x", 1:10),
      neta1 = c(1L, 2L, 2L, 3L, 4L, 4L, 5L, 5L, 6L, 6L),
      neta2 = c(1L, 1L, 2L, 3L, 3L, 4L, 4L, 5L, 5L, 6L),
      fix = FALSE
    )
  )
  .s5 <- c(0L, 0L, 1L, 2L, 1L, 2L)
  expect_true(all(.foceiFdFullSameTie(.p5, rep(FALSE, nrow(.p5)), .s5, paste0("e", 1:6))$free))
  .p6 <- rbind(c(1L, 1L), c(2L, 1L), c(2L, 2L), c(3L, 3L), c(4L, 4L))
  expect_true(all(.foceiFdFullSameTie(.p6, rep(FALSE, 5), c(0L, 0L, 1L, 2L), paste0("e", 1:4))$free))
  # a copy of half a correlated block, and a copy ahead of its master
  .p7 <- rbind(c(1L, 1L), c(2L, 1L), c(2L, 2L), c(3L, 3L))
  expect_true(all(.foceiFdFullSameTie(.p7, rep(FALSE, 4), c(0L, 0L, 1L), .onm)$free))
  expect_true(all(.foceiFdFullSameTie(.p7, rep(FALSE, 4), c(2L, 0L, 0L), .onm)$free))
  # a 2x2 master with its first variance fixed: copies index the free pairs
  .p4 <- rbind(c(1L, 1L), c(2L, 1L), c(2L, 2L), c(3L, 3L), c(4L, 3L), c(4L, 4L))
  .tie <- .foceiFdFullSameTie(.p4, c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), c(0L, 0L, 1L, 2L), paste0("e", 1:4))
  expect_identical(.tie$free, c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_identical(.tie$cpA, c(4L, 4L))
  expect_identical(.tie$cpB, c(3L, 4L))
  expect_identical(.tie$cpOf, c(0L, 1L))
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
  .p <- suppressWarnings(.foceiFdFullParams(.e))
  expect_identical(.p$names, c("tcl", "add.sd", "om.eta.a", "cov.eta.b.eta.a", "om.eta.b"))
  expect_identical(.p$omA, c(1L, 2L, 2L))
  expect_identical(.p$omB, c(1L, 1L, 2L))
  expect_identical(.p$cpA, c(3L, 4L, 4L, 5L, 6L, 6L))
  expect_identical(.p$cpB, c(3L, 3L, 4L, 5L, 5L, 6L))
  expect_identical(.p$cpOf, c(0L, 1L, 2L, 0L, 1L, 2L))
  # every master fixed: no Omega parameter, and the names stay a character vector
  .ui2 <- .ui
  .idf <- .ui2$iniDf
  .idf$fix[!is.na(.idf$neta1) & .idf$neta1 <= 2] <- TRUE
  suppressWarnings(suppressMessages(.ui2$iniDf <- .idf))
  .e$ui <- .ui2
  .p <- suppressWarnings(.foceiFdFullParams(.e))
  expect_identical(.p$names, c("tcl", "add.sd"))
  expect_length(.p$omA, 0L)
})

test_that("the full FD stage moves a same() copy with its master (#1316)", {
  skip_on_cran()
  .f <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      eta.cl ~ 0.3
      eta.v ~ same()
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2(.f, nlmixr2data::theo_sd, "focei", control = foceiControl(print = 0, calcTables = FALSE))
  ))
  .R <- .fit$env$.fdFullR
  expect_identical(rownames(.R), c("tka", "tcl", "tv", "add.sd", "om.eta.cl"))
  # its curvature is the objective's when eta.cl and eta.v move together
  .ofv <- function(v) {
    .ui <- suppressMessages(eval(bquote(rxode2::ini(.fit$finalUi, eta.cl ~ .(v)))))
    suppressWarnings(suppressMessages(nlmixr2(
      .ui,
      nlmixr2data::theo_sd,
      "focei",
      control = foceiControl(
        print = 0,
        calcTables = FALSE,
        maxOuterIterations = 0,
        covMethod = ""
      )
    )))$objf
  }
  .x0 <- .fit$env$.fdFullX0[5]
  .h <- .fit$env$.fdFullH[5]
  .r <- 0.5 * (.ofv(.x0 + .h) - 2 * .ofv(.x0) + .ofv(.x0 - .h)) / .h^2
  expect_equal(.R[5, 5], .r, tolerance = 0.03)
})
