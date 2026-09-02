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

.twoLevelCtl <- function(nBurn = 60, nEm = 80, ...) {
  saemControl(nBurn = nBurn, nEm = nEm, seed = 42L, print = 0L,
              covMethod = "", calcTables = FALSE, ...)
}

test_that("iovMethod='twoLevel' estimates one occasion variance in closed form", {
  skip_on_cran()
  .d <- .twoLevelFitData()
  .f2 <- suppressWarnings(nlmixr2(.twoLevelFitModel(), .d, est = "saem",
                                  control = .twoLevelCtl(iovMethod = "twoLevel")))
  .om <- if (is.list(.f2$omega)) .f2$omega$id else .f2$omega

  # the mechanism: the occasion term is an ordinary omega entry, not a
  # population parameter multiplying a unit-variance eta.  There is no
  # magnitude theta at all.
  expect_false("iov.cl" %in% names(.f2$fixef))
  expect_true(all(c("rx.iov.cl.1", "rx.iov.cl.2") %in% rownames(.om)))

  # the equality constraint (one Psi across occasions) is imposed in the M-step,
  # so the two occasions carry exactly the same variance -- not merely a similar
  # one
  .v1 <- .om["rx.iov.cl.1", "rx.iov.cl.1"]
  .v2 <- .om["rx.iov.cl.2", "rx.iov.cl.2"]
  expect_identical(.v1, .v2)
  expect_true(is.finite(.v1))
  expect_true(.v1 > 0)
})

test_that("the two-level occasion variance does not collapse the way the shared rewrite's does", {
  skip_on_cran()
  .d <- .twoLevelFitData()
  .f0 <- suppressWarnings(nlmixr2(.twoLevelFitModel(), .d, est = "saem",
                                  control = .twoLevelCtl(iovMethod = "theta")))
  .f2 <- suppressWarnings(nlmixr2(.twoLevelFitModel(), .d, est = "saem",
                                  control = .twoLevelCtl(iovMethod = "twoLevel")))
  .om2 <- if (is.list(.f2$omega)) .f2$omega$id else .f2$omega
  .iov0 <- .f0$omega$occ[1, 1]
  .iov2 <- .om2["rx.iov.cl.1", "rx.iov.cl.1"]

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
