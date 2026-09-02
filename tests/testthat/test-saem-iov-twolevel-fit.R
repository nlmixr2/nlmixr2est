# Fit-based tests for saemControl(iovMethod = "twoLevel").  Slow; kept out of
# the essential push/PR subset via .slowBatches in tests/testthat.R.

.twoLevelFitData <- function() {
  .d <- nlmixr2data::theo_md
  .d$occ <- 1
  .d$occ[.d$TIME >= 144] <- 2
  .d
}

# nolint start: object_usage_linter. rxode2 ini()/model() are NSE blocks:
# every assignment here is model specification, consumed by rxode2, not a
# local variable lintr can see used.
.twoLevelFitModel <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    iov.cl ~ 0.1 | occ
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl + iov.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}
# nolint end

.twoLevelCtl <- function(nBurn = 60, nEm = 80, covMethod = "", ...) {
  saemControl(nBurn = nBurn, nEm = nEm, seed = 42L, print = 0L,
              covMethod = covMethod, calcTables = FALSE, ...)
}

test_that("iovMethod='twoLevel' estimates one occasion variance in closed form", {
  skip_on_cran()
  .d <- .twoLevelFitData()
  .f2 <- suppressWarnings(nlmixr2(.twoLevelFitModel(), .d, est = "saem",
                                  control = .twoLevelCtl(iovMethod = "twoLevel")))

  # the mechanism: the occasion term is an ordinary omega entry, not a
  # population parameter multiplying a unit-variance eta.  There is no
  # magnitude theta at all.
  expect_false("iov.cl" %in% names(.f2$fixef))

  # and the fit presents the same way the shared rewrite's does: $omega split
  # into $id and $occ, the user's `iov.cl ~ v | occ` row restored, an $iov
  # deviation table, and no rx.* leaking into the fit
  expect_true(is.list(.f2$omega))
  expect_equal(names(.f2$omega), c("id", "occ"))
  expect_equal(rownames(.f2$omega$occ), "iov.cl")
  expect_true(is.finite(.f2$omega$occ[1, 1]))
  expect_true(.f2$omega$occ[1, 1] > 0)
  .ini <- .f2$ui$iniDf
  expect_equal(.ini$condition[.ini$name == "iov.cl"], "occ")
  expect_false(any(grepl("^rx[.]", .ini$name)))
  expect_true(all(c("ID", "occ", "iov.cl") %in% names(.f2$iov$occ)))
  expect_length(grep("^rx[.]", names(.f2)), 0L)
})

test_that("the two-level occasion variance does not collapse the way the shared rewrite's does", {
  skip_on_cran()
  .d <- .twoLevelFitData()
  .f0 <- suppressWarnings(nlmixr2(.twoLevelFitModel(), .d, est = "saem",
                                  control = .twoLevelCtl(iovMethod = "theta")))
  .f2 <- suppressWarnings(nlmixr2(.twoLevelFitModel(), .d, est = "saem",
                                  control = .twoLevelCtl(iovMethod = "twoLevel")))
  .iov0 <- .f0$omega$occ[1, 1]
  .iov2 <- .f2$omega$occ[1, 1]

  # the shared rewrite estimates this variance through the annealed phi0 path
  # and drives it toward zero; the closed-form M-step does not
  expect_true(.iov2 > 2 * .iov0)

  # the structural parameters are unaffected -- this is the same model, fitted
  # better, not a different one
  expect_equal(.f2$fixef[["tka"]], .f0$fixef[["tka"]], tolerance = 0.1)
  expect_equal(.f2$fixef[["tcl"]], .f0$fixef[["tcl"]], tolerance = 0.02)
  expect_equal(.f2$fixef[["tv"]], .f0$fixef[["tv"]], tolerance = 0.01)
  expect_equal(.f2$fixef[["add.sd"]], .f0$fixef[["add.sd"]], tolerance = 0.05)
})

test_that("a model outside the two-level scope falls back to the shared rewrite", {
  skip_on_cran()
  .d <- .twoLevelFitData()
  .d$occ2 <- ifelse(.d$TIME %% 2 < 1, 1, 2)
  # nolint start: object_usage_linter. rxode2 ini()/model() are NSE blocks:
  # every assignment here is model specification, consumed by rxode2, not a
  # local variable lintr can see used.
  .two <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      iov.cl ~ 0.1 | occ
      iov.v ~ 0.1 | occ2
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv + eta.v + iov.v)
      linCmt() ~ add(add.sd)
    })
  }
  # nolint end
  .f <- suppressWarnings(nlmixr2(.two(), .d, est = "saem",
                                 control = .twoLevelCtl(nBurn = 15, nEm = 15,
                                                        iovMethod = "twoLevel")))
  # it fitted, by the shared rewrite, and said so
  expect_true(inherits(.f, "nlmixr2FitCore"))
  expect_true(any(grepl("two-level IOV needs one occasion variable", .f$runInfo)))
  expect_true(is.list(.f$omega))
  expect_true(all(c("occ", "occ2") %in% names(.f$omega)))
})

test_that("a two-level fit's covariance carries one row for the occasion variance", {
  skip_on_cran()
  .d <- .twoLevelFitData()
  .f <- suppressWarnings(nlmixr2(.twoLevelFitModel(), .d, est = "saem",
                                 control = .twoLevelCtl(covMethod = "linFim",
                                                        iovMethod = "twoLevel")))
  .cv <- .f$cov
  skip_if(is.null(.cv) || !is.matrix(.cv))
  # the K per-occasion columns are ONE parameter; they must not appear
  # separately, and certainly not under their internal names
  expect_length(grep("^om[.]rx[.]", rownames(.cv)), 0L)
  expect_true("om.iov.cl" %in% rownames(.cv))
  expect_equal(rownames(.cv), colnames(.cv))
  expect_true(is.finite(.cv["om.iov.cl", "om.iov.cl"]))
  expect_true(.cv["om.iov.cl", "om.iov.cl"] > 0)
})

test_that("iovMethod='collapsed' imposes both constraints exactly", {
  skip_on_cran()
  .d <- .twoLevelFitData()
  .f <- suppressWarnings(nlmixr2(.twoLevelFitModel(), .d, est = "saem",
                                 control = .twoLevelCtl(iovMethod = "collapsed")))

  # the fit presents in the user's own parameterization, like every other path
  expect_true("tcl" %in% names(.f$fixef))
  expect_false(any(grepl("^rx[.]", names(.f$fixef))))
  expect_equal(names(.f$omega), c("id", "occ"))
  expect_true("eta.cl" %in% rownames(.f$omega$id))
  expect_equal(rownames(.f$omega$occ), "iov.cl")
  expect_length(grep("^rx[.]", names(.f)), 0L)
  expect_length(grep("^rx[.]", names(.f$env$ranef)), 0L)
  expect_true("eta.cl" %in% names(.f$env$ranef))

  # for K = 2 the deviations are exactly antisymmetric within a subject,
  # because b_i is defined as the subject's mean over occasions
  .iov <- .f$iov[[1]]
  expect_true(all(c("ID", "occ", "iov.cl") %in% names(.iov)))
  .s <- tapply(.iov$iov.cl, .iov$ID, sum)
  expect_true(max(abs(.s)) < 1e-8)

  # both variances are positive and finite
  expect_true(is.finite(.f$omega$occ[1, 1]) && .f$omega$occ[1, 1] > 0)
  expect_true(.f$omega$id["eta.cl", "eta.cl"] > 0)
})

test_that("the collapsed and two-level paths do not contaminate each other", {
  skip_on_cran()
  # .uiIovEnv is process-global, and each path installs its own post-fit
  # restoration; a stale marker sends a fit through the wrong one
  .d <- .twoLevelFitData()
  .fc <- suppressWarnings(nlmixr2(.twoLevelFitModel(), .d, est = "saem",
                                  control = .twoLevelCtl(nBurn = 15, nEm = 15,
                                                         iovMethod = "collapsed")))
  .ft <- suppressWarnings(nlmixr2(.twoLevelFitModel(), .d, est = "saem",
                                  control = .twoLevelCtl(nBurn = 15, nEm = 15,
                                                         iovMethod = "twoLevel")))
  .fl <- suppressWarnings(nlmixr2(.twoLevelFitModel(), .d, est = "saem",
                                  control = .twoLevelCtl(nBurn = 15, nEm = 15,
                                                         iovMethod = "theta")))
  for (.f in list(.fc, .ft, .fl)) {
    expect_equal(names(.f$omega), c("id", "occ"))
    expect_equal(rownames(.f$omega$occ), "iov.cl")
    expect_true("tcl" %in% names(.f$fixef))
    expect_length(grep("^rx[.]", names(.f)), 0L)
  }
})
