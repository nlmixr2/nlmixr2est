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

test_that("fix() on an occasion parameter is respected under 'omega'", {
  skip_on_cran()
  # `iov.cl ~ fix(0.1) | occ` fixes the VARIANCE.  Under "theta" that
  # rides on the magnitude theta; under "omega" the variance IS the omega
  # block, so the flag has to land there -- otherwise the parameter is
  # quietly estimated while the fit still reports it as fixed.
  .f <- .sameMod("iov.cl ~ fix(0.1) | occ",
                 "cl <- exp(tcl + iov.cl)\n v <- exp(tv)")
  .ini <- .applyIov(.f, "omega")$ui$iniDf
  .occ <- .ini[grepl("^rx\\.iov\\.cl\\.", .ini$name), ]
  expect_true(all(.occ$fix))
  expect_equal(.occ$est, rep(0.1, 2))

  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.f(), .sameData(), "focei",
                        foceiControl(print = 0, covMethod = "",
                                     iovMethod = "omega"))))
  .row <- .fit$ui$iniDf[.fit$ui$iniDf$name == "iov.cl", ]
  expect_true(.row$fix)
  expect_equal(.row$est, 0.1)
})

test_that("a fixed correlated block stays a repeated block", {
  skip_on_cran()
  # every occasion carries the same `fix` flags, so lotri still re-emits
  # `same()` -- a mismatch there would silently dissolve the repetition
  .f <- .sameMod("iov.cl + iov.v ~ fix(0.1,\n 0.03, 0.2) | occ",
                 "cl <- exp(tcl + iov.cl)\n v <- exp(tv + iov.v)")
  .ini <- .applyIov(.f, "omega")$ui$iniDf
  .occ <- .ini[grepl("^\\(?rx\\.iov\\.", .ini$name), ]
  expect_true(all(.occ$fix))
  expect_true(any(grepl(":same:", .ini$condition, fixed = TRUE)))
})

test_that("the two expansions agree exactly where they coincide", {
  skip_on_cran()
  # At an occasion variance of 1 the parameterizations are literally the
  # same problem -- "theta" scales a unit-variance eta by a theta of 1,
  # "omega" gives the eta a variance of 1 -- so the objective must agree
  # to machine precision.  This is the claim NEWS.md makes; if it ever
  # stops holding, the expansions have diverged.
  .f <- .sameMod("iov.cl ~ 1 | occ", "cl <- exp(tcl + iov.cl)\n v <- exp(tv)")
  .obj <- vapply(c("theta", "omega"), function(.m) {
    suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(.f(), .sameData(), "focei",
                          foceiControl(print = 0, sigdig = 7, covMethod = "",
                                       maxOuterIterations = 0L,
                                       maxInnerIterations = 100000L,
                                       iovMethod = .m))))$objf
  }, double(1))
  expect_equal(.obj[[1]], .obj[[2]], tolerance = 1e-10)
})

test_that("analytic covariance bows out for an 'omega' IOV fit", {
  skip_on_cran()
  # The analytic IOV branch reads the magnitude theta as the occasion
  # standard deviation.  That theta is FIXED AT ONE here, so it would
  # overwrite the estimated per-occasion variances rather than merely be
  # conservative.  The mode is detected from the occasion etas -- the
  # control still says "auto" after resolving, and with only ONE occasion
  # there is no repeated block for `omegaSameMap` to report either.
  .d <- .sameData()
  .d$occ <- 1L
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.corMod(), .d, "focei",
                        foceiControl(print = 0, covMethod = "analytic"))))
  expect_false(identical(.fit$covMethod, "analytic"))
  # the occasion variances are the ESTIMATED ones, not the theta's 1
  expect_true(all(diag(.fit$omega$occ) != 1))
})

test_that("only a method that honours the shared block may use it", {
  skip_on_cran()
  expect_true(.isIovSameMethod("focei"))
  expect_true(.isIovSameMethod("laplace"))
  # SAEM's omega is a moment M-step with no notion of a repeated block;
  # the variational and nonparametric methods build their own omega
  expect_false(.isIovSameMethod("saem"))
  expect_false(.isIovSameMethod("npag"))
  expect_false(.isIovSameMethod("fbvi"))
  expect_false(.isIovSameMethod(NA_character_))
  expect_false(.isIovSameMethod(character(0)))
})

test_that("saem still refuses a correlated occasion block", {
  skip_on_cran()
  # This is the one that would go WRONG rather than merely unsupported:
  # SAEM would estimate each occasion's block independently and
  # `.uiFinalizeIov()` would report occasion one and discard the rest.
  # `assertRxUiIovNoCor()` cannot catch it -- by the time SAEM runs, the
  # block has already been expanded to `id`/`id:same:` rows whose BASE
  # condition is "id" -- so `"auto"` must never choose "omega" here.
  expect_error(
    .uiApplyIov(rxode2::rxode2(.corMod()), "saem", .sameData(),
                saemControl(print = 0)),
    "not supported by 'saem'")
  # asking for it outright is refused by name, not silently downgraded
  expect_error(
    .uiApplyIov(rxode2::rxode2(.corMod()), "saem", .sameData(),
                foceiControl(iovMethod = "omega")),
    "does not honour")
})

test_that("the FOCEi family really shares the repeated block", {
  skip_on_cran()
  # Running is not enough: the reported `$occ` is read off occasion one,
  # so two independently estimated blocks would look fine.  The
  # parameter count is the tell -- 4 thetas + eta.ka + ONE 2x2 (3) = 8,
  # against 11 if each occasion were estimated on its own.
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.corMod(), .sameData(), "focei",
                        foceiControl(print = 0, covMethod = ""))))
  expect_equal(attr(logLik(.fit), "df"), 8L)
})
