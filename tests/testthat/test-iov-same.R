# `iovMethod = "omega"`: the occasion blocks are one estimated covariance
# repeated per occasion (NONMEM's `$OMEGA BLOCK(n) SAME`), which is what
# lets inter-occasion parameters be CORRELATED.

.sameData <- function() {
  .d <- nlmixr2data::theo_sd
  .d$occ <- 1 + (.d$TIME >= 5)
  .d
}

.sameMod <- function(blk, mdl) {
  eval(parse(text = sprintf('function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- c(0, 0.7)
          eta.ka ~ 0.6
          %s })
    model({ ka <- exp(tka + eta.ka)
            %s
            linCmt() ~ add(add.sd) }) }', blk, mdl)))
}

.corMod <- function() {
  .sameMod("iov.cl + iov.v ~ c(0.1,\n 0.03, 0.2) | occ",
           "cl <- exp(tcl + iov.cl)\n v <- exp(tv + iov.v)")
}

.diagMod <- function() {
  .sameMod("iov.cl ~ 0.1 | occ",
           "cl <- exp(tcl + iov.cl)\n v <- exp(tv)")
}

.applyIov <- function(f, m) {
  .uiApplyIov(rxode2::rxode2(f()), "focei", .sameData(),
              foceiControl(iovMethod = m))
}

test_that("'auto' picks the expansion from whether the block is correlated", {
  skip_on_cran()
  # a diagonal block keeps the long standing "theta" expansion: unit
  # variance etas scaled by a magnitude theta
  .ini <- .applyIov(.diagMod, "auto")$ui$iniDf
  expect_false(any(grepl(":same:", .ini$condition, fixed = TRUE)))
  expect_equal(.ini$est[.ini$name == "iov.cl"], sqrt(0.1))
  expect_true(all(.ini$fix[grepl("^rx\\.iov\\.cl\\.", .ini$name)]))
  # a correlated block cannot be represented that way, so it routes to the
  # shared omega block instead
  .ini <- .applyIov(.corMod, "auto")$ui$iniDf
  expect_true(any(grepl(":same:", .ini$condition, fixed = TRUE)))
})

test_that("iovMethod='omega' expands occasion-major into repeated blocks", {
  skip_on_cran()
  .ui <- .applyIov(.corMod, "omega")$ui
  .ini <- .ui$iniDf
  .eta <- .ini[!is.na(.ini$neta1), ]

  # the magnitude theta is kept, but fixed at one -- the variability is in
  # the omega block now
  .th <- .ini[!is.na(.ini$ntheta), ]
  expect_equal(.th$est[.th$name == "iov.cl"], 1)
  expect_equal(.th$est[.th$name == "iov.v"], 1)
  expect_true(all(.th$fix[.th$name %in% c("iov.cl", "iov.v")]))

  # OCCASION-major: each occasion's block is contiguous, so a repeated
  # block sits immediately after the block it repeats
  .diag <- .eta[.eta$neta1 == .eta$neta2, ]
  expect_equal(.diag$name,
               c("eta.ka", "rx.iov.cl.1", "rx.iov.v.1",
                 "rx.iov.cl.2", "rx.iov.v.2"))

  # occasion one IS the block; the rest point back at it by eta NAME
  expect_equal(.eta$condition[.eta$name == "rx.iov.cl.1"], "id")
  expect_equal(.eta$condition[.eta$name == "rx.iov.cl.2"],
               "id:same:rx.iov.cl.1")
  expect_equal(.eta$condition[.eta$name == "rx.iov.v.2"],
               "id:same:rx.iov.v.1")
  # ... covariances included, which is the whole point
  expect_equal(.eta$condition[.eta$name == "(rx.iov.cl.2,rx.iov.v.2)"],
               "id:same:rx.iov.cl.1:rx.iov.v.1")
  expect_equal(.eta$est[.eta$name == "(rx.iov.cl.2,rx.iov.v.2)"], 0.03)

  # A per-occasion label would differ between a block and its copy, and
  # lotri only re-emits `same()` when values, `fix` AND labels match -- so
  # a label here silently breaks the repetition on the next `$fun()`.
  expect_true(all(is.na(.eta$label[grepl("^rx\\.iov\\.", .eta$name)])))
})

test_that("the expanded omega is one matrix of identical blocks", {
  skip_on_cran()
  .ui <- .applyIov(.corMod, "omega")$ui
  .om <- .ui$omega
  expect_equal(dim(.om), c(5L, 5L))
  # occasion 1 block == occasion 2 block, correlation and all
  expect_equal(unname(.om[2:3, 2:3]), unname(.om[4:5, 4:5]))
  expect_equal(unname(.om[2:3, 2:3]),
               matrix(c(0.1, 0.03, 0.03, 0.2), 2, 2))
  # ... and the repetition is recorded, by eta index
  expect_equal(.ui$omegaSameMap, c(0L, 0L, 0L, 2L, 3L))
})

test_that("a repeated block shares its master's cholesky parameters", {
  skip_on_cran()
  .ui <- .applyIov(.corMod, "omega")$ui
  .om <- .ui$omega
  .free <- rxode2::rxSymInvCholCreate(.om, "sqrt", same = .ui$omegaSameMap)
  .all <- rxode2::rxSymInvCholCreate(.om, "sqrt")
  # eta.ka (1) + ONE 2x2 block (3) rather than eta.ka + TWO blocks (6)
  expect_equal(length(.free$theta), 4L)
  expect_equal(length(.all$theta), 7L)
})

test_that("iovMethod='omega' is allowed on a diagonal block too", {
  skip_on_cran()
  # not just permitted but required: it is what makes the two expansions
  # numerically comparable on the same model
  .ini <- .applyIov(.diagMod, "omega")$ui$iniDf
  expect_true(any(grepl(":same:", .ini$condition, fixed = TRUE)))
  expect_equal(.ini$est[.ini$name == "iov.cl"], 1)
  expect_equal(.ini$est[.ini$name == "rx.iov.cl.1"], 0.1)
  expect_equal(.ini$est[.ini$name == "rx.iov.cl.2"], 0.1)
})

test_that("iovMethod='theta' still refuses a correlated block", {
  skip_on_cran()
  expect_error(.applyIov(.corMod, "theta"), "correlated inter-occasion")
})

test_that("a correlated IOV model fits and restores the user's own block", {
  skip_on_cran()
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.corMod(), .sameData(), "focei",
                        foceiControl(print = 0, covMethod = ""))))
  .ini <- .fit$ui$iniDf

  # the fit reports the model the USER wrote, not the expansion: one
  # occasion block, off diagonal included, and no `rx.` etas left over
  .eta <- .ini[!is.na(.ini$neta1), ]
  expect_equal(sort(.eta$name[.eta$condition == "occ"]),
               sort(c("iov.cl", "iov.v", "(iov.cl,iov.v)")))
  expect_false(any(grepl("^rx\\.", .ini$name)))
  # the magnitude thetas are gone from the reported parameters
  expect_false(any(c("iov.cl", "iov.v") %in%
                     .ini$name[!is.na(.ini$ntheta)]))

  # the estimated occasion block really is correlated
  .om <- .fit$omega$occ
  expect_equal(dim(.om), c(2L, 2L))
  expect_true(abs(.om[1, 2]) > 0)
  expect_equal(.om[1, 2], .om[2, 1])

  # the reported CV comes from the ESTIMATED variance, not from the
  # magnitude theta -- which is fixed at one, so it would report a
  # constant 131% for every model
  .bck <- grep("Back", names(.fit$parFixedDf))
  expect_equal(.fit$parFixedDf["iov.cl", .bck],
               nlmixr2iovSdCv(sqrt(.om[1, 1])))
  expect_equal(.fit$parFixedDf["iov.v", .bck],
               nlmixr2iovSdCv(sqrt(.om[2, 2])))

  # per-occasion etas are reported against a real occasion, and are NOT
  # rescaled by the magnitude (it is one)
  expect_true(all(!is.na(.fit$iov$occ$occ)))
  expect_equal(sort(unique(.fit$iov$occ$occ)), c(1L, 2L))
})

test_that("two occasion parameters get their occasion number, not NA", {
  skip_on_cran()
  # `.dt`'s occasion column is built from the FIRST occasion parameter's
  # `rx.<name>.<occ>` spellings; stripping the LAST one's prefix left every
  # occasion NA (and warned "NAs introduced by coercion") as soon as a
  # level carried two parameters.  Reproduces on the "theta" path too.
  .f <- .sameMod("iov.cl ~ 0.1 | occ; iov.v ~ 0.2 | occ",
                 "cl <- exp(tcl + iov.cl)\n v <- exp(tv + iov.v)")
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.f(), .sameData(), "focei",
                        foceiControl(print = 0, covMethod = "",
                                     iovMethod = "theta"))))
  expect_true(all(!is.na(.fit$iov$occ$occ)))
  expect_equal(sort(unique(.fit$iov$occ$occ)), c(1L, 2L))
  expect_false(any(grepl("NAs introduced", .fit$runInfo)))
})
