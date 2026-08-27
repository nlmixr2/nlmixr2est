# .uiApplyIov() builds the IOV magnitude theta and the per-occasion etas by
# COPYING an existing iniDf row as a template.  Everything the template
# carries that is per-parameter has to be overwritten, or it silently
# becomes this parameter's -- `prior` was not, and `fix` was read from a
# vector over the whole occasion rather than from the parameter's own row.
#
# These exercise the rewrite directly (no fit), so they stay fast.

.iovData <- function() {
  .d <- nlmixr2data::theo_md
  .d$occ <- 1
  .d$occ[.d$TIME >= 144] <- 2
  .d
}

.iovApply <- function(f, data = .iovData()) {
  .uiApplyIov(rxode2::rxode2(f), "focei", data, foceiControl())
}

test_that("the magnitude theta carries prior(iov.x), not theta #1's", {
  skip_on_cran()
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      prior(tka) ~ dnorm(-7.25, 0.125)   # deliberately distinctive
      prior(iov.cl) ~ dcauchy(0, 1)
      eta.cl ~ 0.3
      iov.cl ~ 0.1 | occ
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  skip_if_not("prior" %in% names(rxode2::rxode2(.mod)$iniDf),
              "this lotri has no prior support")
  .ini <- .iovApply(.mod)$ui$iniDf
  .mag <- .ini[!is.na(.ini$ntheta) & .ini$name == "iov.cl", ]
  expect_equal(nrow(.mag), 1L)
  # the declared prior survives the rewrite that deletes the row it was
  # written on ...
  expect_equal(.mag$prior, "dcauchy(0, 1)")
  # ... and the first theta's prior stays on the first theta
  expect_equal(.ini$prior[.ini$name == "tka"], "dnorm(-7.25, 0.125)")
  # the per-occasion etas are not given a prior at all
  .occ <- .ini[grepl("^rx\\.iov\\.cl\\.", .ini$name), ]
  expect_equal(nrow(.occ), 2L)
  expect_true(all(is.na(.occ$prior)))
})

test_that("an undeclared magnitude, and a fixed one, carry no prior", {
  skip_on_cran()
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      prior(tka) ~ dnorm(-7.25, 0.125)
      eta.cl ~ 0.3
      iov.cl ~ 0.1 | occ
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  skip_if_not("prior" %in% names(rxode2::rxode2(.mod)$iniDf),
              "this lotri has no prior support")
  .ini <- .iovApply(.mod)$ui$iniDf
  expect_true(is.na(.ini$prior[!is.na(.ini$ntheta) & .ini$name == "iov.cl"]))

  # a FIXED magnitude is a constant: it cannot hold a prior, and inheriting
  # the template's made rxode2 refuse the rewritten model outright
  .fixMod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      prior(tka) ~ dnorm(-7.25, 0.125)
      eta.cl ~ 0.3
      iov.cl ~ fix(0.1) | occ
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  .iniF <- .iovApply(.fixMod)$ui$iniDf
  .magF <- .iniF[!is.na(.iniF$ntheta) & .iniF$name == "iov.cl", ]
  expect_true(.magF$fix)
  expect_true(is.na(.magF$prior))
})

test_that("a prior on the FIRST eta does not become the occasion etas'", {
  skip_on_cran()
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      prior(eta.cl) ~ invWishart(4)
      eta.cl ~ 0.3
      iov.cl ~ 0.1 | occ
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  skip_if_not("prior" %in% names(rxode2::rxode2(.mod)$iniDf),
              "this lotri has no prior support")
  .ini <- .iovApply(.mod)$ui$iniDf
  # eta.cl keeps its own prior; the fixed per-occasion etas copied from it
  # get none (rxode2 refuses a prior on a fixed parameter, so inheriting it
  # made the whole model unusable)
  expect_equal(.ini$prior[.ini$name == "eta.cl"], "invWishart(4)")
  expect_true(all(is.na(.ini$prior[grepl("^rx\\.iov\\.cl\\.", .ini$name)])))
})

test_that("several occasion parameters ride one occasion variable", {
  skip_on_cran()
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      eta.cl ~ 0.3
      iov.cl ~ 0.1 | occ
      iov.v ~ fix(0.04) | occ
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv + iov.v)
      linCmt() ~ add(add.sd)
    })
  }
  .ini <- .iovApply(.mod)$ui$iniDf
  .th <- .ini[!is.na(.ini$ntheta), ]
  # one magnitude theta each, ONCE (the occasion variable used to be
  # visited once per parameter riding it, duplicating every theta) ...
  expect_equal(sum(.th$name == "iov.cl"), 1L)
  expect_equal(sum(.th$name == "iov.v"), 1L)
  # ... and `fix` comes from each parameter's own row rather than from a
  # vector over the whole occasion, which errored with "replacement has 2
  # rows, data has 1"
  expect_false(.th$fix[.th$name == "iov.cl"])
  expect_true(.th$fix[.th$name == "iov.v"])
  expect_equal(.th$est[.th$name == "iov.cl"], sqrt(0.1))
  expect_equal(.th$est[.th$name == "iov.v"], sqrt(0.04))
  # per-occasion etas for both, two occasions each
  expect_equal(sum(grepl("^rx\\.iov\\.cl\\.", .ini$name)), 2L)
  expect_equal(sum(grepl("^rx\\.iov\\.v\\.", .ini$name)), 2L)
})

test_that("correlated occasion random effects are refused explicitly", {
  skip_on_cran()
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      eta.cl ~ 0.3
      iov.cl + iov.v ~ c(0.1,
                         0.01, 0.1) | occ
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv + iov.v)
      linCmt() ~ add(add.sd)
    })
  }
  # the expansion gives each occasion parameter its own magnitude theta and
  # unit-variance etas, which cannot carry a correlation; the off-diagonal
  # row used to be treated as a third occasion parameter, giving a syntax
  # error from rxRename() about "rx.(iov.cl,iov.v)="
  expect_error(.iovApply(.mod), "correlated inter-occasion")
})

test_that("an occasion parameter declared twice is named in the error", {
  skip_on_cran()
  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
      eta.cl ~ 0.3
      iov.cl ~ 0.1 | occ
      iov.cl ~ 0.15 | occ
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  # rxode2 builds this ui (two eta rows share the name), so the rewrite is
  # what has to catch it; every per-parameter field would otherwise be a
  # vector and it would die on "replacement has 2 rows, data has 1" without
  # saying which parameter
  expect_equal(sum(rxode2::rxode2(.mod)$iniDf$name == "iov.cl"), 2L)
  expect_error(.iovApply(.mod), "'iov.cl' has 2 variance declarations")
})
