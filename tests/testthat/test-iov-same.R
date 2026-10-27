# `iovMethod = "omega"`: the occasion blocks are one estimated covariance
# repeated per occasion (NONMEM's `$OMEGA BLOCK(n) SAME`), which is what
# lets inter-occasion parameters be CORRELATED.

.sameData <- function() {
  .d <- nlmixr2data::theo_sd
  .d$occ <- 1 + (.d$TIME >= 5)
  .d
}

.sameMod <- function(blk, mdl) {
  eval(parse(
    text = sprintf(
      'function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- c(0, 0.7)
          eta.ka ~ 0.6
          %s })
    model({ ka <- exp(tka + eta.ka)
            %s
            linCmt() ~ add(add.sd) }) }',
      blk,
      mdl
    )
  ))
}

.corMod <- function() {
  .sameMod("iov.cl + iov.v ~ c(0.1,\n 0.03, 0.2) | occ", "cl <- exp(tcl + iov.cl)\n v <- exp(tv + iov.v)")
}

.diagMod <- function() {
  .sameMod("iov.cl ~ 0.1 | occ", "cl <- exp(tcl + iov.cl)\n v <- exp(tv)")
}

.applyIov <- function(f, m) {
  .uiApplyIov(rxode2::rxode2(f()), "focei", .sameData(), foceiControl(iovMethod = m))
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
  expect_equal(.diag$name, c("eta.ka", "rx.iov.cl.1", "rx.iov.v.1", "rx.iov.cl.2", "rx.iov.v.2"))

  # occasion one IS the block; the rest point back at it by eta NAME
  expect_equal(.eta$condition[.eta$name == "rx.iov.cl.1"], "id")
  expect_equal(.eta$condition[.eta$name == "rx.iov.cl.2"], "id:same:rx.iov.cl.1")
  expect_equal(.eta$condition[.eta$name == "rx.iov.v.2"], "id:same:rx.iov.v.1")
  # ... covariances included, which is the whole point
  expect_equal(.eta$condition[.eta$name == "(rx.iov.cl.2,rx.iov.v.2)"], "id:same:rx.iov.cl.1:rx.iov.v.1")
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
  expect_equal(unname(.om[2:3, 2:3]), matrix(c(0.1, 0.03, 0.03, 0.2), 2, 2))
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
    nlmixr2est::nlmixr2(.corMod(), .sameData(), "focei", foceiControl(print = 0, covMethod = ""))
  ))
  .ini <- .fit$ui$iniDf

  # the fit reports the model the USER wrote, not the expansion: one
  # occasion block, off diagonal included, and no `rx.` etas left over
  .eta <- .ini[!is.na(.ini$neta1), ]
  expect_equal(sort(.eta$name[.eta$condition == "occ"]), sort(c("iov.cl", "iov.v", "(iov.cl,iov.v)")))
  expect_false(any(grepl("^rx\\.", .ini$name)))
  # the magnitude thetas are gone from the reported parameters
  expect_false(any(
    c("iov.cl", "iov.v") %in%
      .ini$name[!is.na(.ini$ntheta)]
  ))

  # the estimated occasion block really is correlated
  .om <- .fit$omega$occ
  expect_equal(dim(.om), c(2L, 2L))
  expect_true(abs(.om[1, 2]) > 0)
  expect_equal(.om[1, 2], .om[2, 1])

  # the reported CV comes from the ESTIMATED variance, not from the
  # magnitude theta -- which is fixed at one, so it would report a
  # constant 131% for every model
  .bck <- grep("Back", names(.fit$parFixedDf))
  expect_equal(.fit$parFixedDf["iov.cl", .bck], nlmixr2iovSdCv(sqrt(.om[1, 1])))
  expect_equal(.fit$parFixedDf["iov.v", .bck], nlmixr2iovSdCv(sqrt(.om[2, 2])))

  # per-occasion etas are reported against a real occasion, and are NOT
  # rescaled by the magnitude (it is one)
  expect_true(all(!is.na(.fit$iov$occ$occ)))
  expect_equal(sort(unique(.fit$iov$occ$occ)), c(1L, 2L))
})

test_that("an id-level same() block beside IOV fits in both expansions", {
  skip_on_cran()
  # eta.v repeats eta.cl's block at the id level (condition "id:same:eta.cl"),
  # so it is a between-subject eta; iov.cl is the only occasion parameter
  .f <- .sameMod(
    "eta.cl ~ 0.3; eta.v ~ same(); iov.cl ~ 0.1 | occ",
    "cl <- exp(tcl + eta.cl + iov.cl)\n v <- exp(tv + eta.v)"
  )
  for (.m in c("theta", "omega")) {
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(.f(), .sameData(), "focei", foceiControl(print = 0, covMethod = "", iovMethod = .m))
    ))
    .ini <- .fit$ui$iniDf
    .eta <- .ini[!is.na(.ini$neta1), ]
    expect_identical(.eta$name, c("eta.ka", "eta.cl", "eta.v", "iov.cl"), label = .m)
    expect_identical(.eta$condition[.eta$name == "iov.cl"], "occ", label = .m)
    expect_false(any(grepl("^rx\\.", .ini$name)), label = .m)
    expect_false("iov.cl" %in% .ini$name[!is.na(.ini$ntheta)], label = .m)
    # the copy kept its master's variance
    expect_identical(dimnames(.fit$omega$id)[[1]], c("eta.ka", "eta.cl", "eta.v"), label = .m)
    expect_equal(.fit$omega$id["eta.v", "eta.v"], .fit$omega$id["eta.cl", "eta.cl"], label = .m)
    expect_identical(dimnames(.fit$omega$occ)[[1]], "iov.cl", label = .m)
    .bck <- grep("Back", names(.fit$parFixedDf))
    expect_equal(.fit$parFixedDf["iov.cl", .bck], nlmixr2iovSdCv(sqrt(.fit$omega$occ[1, 1])), label = .m)
    expect_equal(sort(unique(.fit$iov$occ$occ)), c(1L, 2L), label = .m)
    expect_false(anyNA(.fit$iov$occ$iov.cl), label = .m)
  }
})

# an id-level same() block and no IOV: eta.v repeats eta.cl's block
.idSameMod <- function() {
  .sameMod("eta.cl ~ 0.3; eta.v ~ same()", "cl <- exp(tcl + eta.cl)\n v <- exp(tv + eta.v)")
}

.sameRefusal <- function(est) {
  paste0(
    "'",
    est,
    "' does not honour a repeated (same()) omega block\n",
    "use fo, foi, foce, focei, focep, laplace, agq, emvi, fbvi or vae"
  )
}

# register a do-nothing nlmixr2Est method `name` in `ns`, carrying
# attr(, attrName) <- value
.fakeEst <- function(env, ...) NULL
.regFake <- function(ns, name, attrName, value) {
  .fn <- .fakeEst
  attr(.fn, attrName) <- value
  registerS3method("nlmixr2Est", name, .fn, envir = ns)
}

.errMsg <- function(expr) {
  tryCatch(
    {
      expr
      NULL
    },
    error = conditionMessage
  )
}

test_that("every method that would drop a same() tie refuses the model", {
  skip_on_cran()
  .ui <- rxode2::rxode2(.idSameMod())
  # nlmixr2est's own methods; a method another test registers is not
  # nlmixr2est's to classify
  .ests <- Filter(
    function(e) exists(paste0("nlmixr2Est.", e), envir = asNamespace("nlmixr2est"), inherits = FALSE),
    nlmixr2AllEst()
  )
  .refused <- .ests[!vapply(.ests, .isSameBlockMethod, logical(1), USE.NAMES = FALSE)]
  # measured on a same() model: nlme drops the tie; saem and the imp, qrpem and
  # npag families crash on it; the optimizers refuse any eta
  expect_identical(
    sort(.refused),
    sort(c(
      "nlme",
      "saem",
      "imp",
      "impmap",
      "qrpem",
      "npag",
      "npb",
      "inpag",
      "inpb",
      "mnpag",
      "mnpb",
      "bfgs",
      "bobyqa",
      "brent",
      "cg",
      "lbfgsb",
      "lbfgsb3c",
      "n1qn1",
      "neldermead",
      "newuoa",
      "nlm",
      "nlminb",
      "nls",
      "optim",
      "sann",
      "trust",
      "uobyqa"
    ))
  )
  for (.e in .refused) {
    expect_identical(.errMsg(.nlmixr0preProcessAssertSameBlock(.ui, .e, NULL, NULL)), .sameRefusal(.e))
  }
  for (.e in setdiff(.ests, .refused)) {
    expect_null(.nlmixr0preProcessAssertSameBlock(.ui, .e, NULL, NULL), label = .e)
  }
  # an occasion block without same() is no repeated block
  .iov <- rxode2::rxode2(.diagMod())
  for (.e in .ests) {
    expect_null(.nlmixr0preProcessAssertSameBlock(.iov, .e, NULL, NULL), label = .e)
  }
})

test_that("the refusal asks nothing of an est with no method, or an eta row with no condition", {
  skip_on_cran()
  # foceiLikLoad() runs the hooks with whatever est it is given (the load
  # itself is in the weekly test-iov-same-fit.R); with no nlmixr2Est method
  # there is nothing to refuse
  expect_null(.nlmixr0preProcessAssertSameBlock(rxode2::rxode2(.idSameMod()), "x", NULL, NULL))
  # an eta row whose condition is missing is no repeated block
  .ui <- list(iniDf = data.frame(neta1 = c(1, 2), condition = c("id", NA_character_)))
  expect_null(.nlmixr0preProcessAssertSameBlock(.ui, "saem", NULL, NULL))
})

test_that("the refusal runs before the IOV hooks and reaches the user", {
  skip_on_cran()
  expect_identical(preProcessHooks()[1], ".nlmixr0preProcessAssertSameBlock")
  expect_identical(
    .errMsg(suppressMessages(nlmixr2est::nlmixr2(.idSameMod(), .sameData(), "saem", saemControl(print = 0)))),
    .sameRefusal("saem")
  )
})

test_that("an out-of-tree method keeps a same() block through its attribute", {
  skip_on_cran()
  .ns <- asNamespace("nlmixr2est")
  withr::defer(rm(
    list = paste0("nlmixr2Est.", c("fakeSame", "fakeIovSame", "fakeFixed", "fakeNone")),
    envir = .ns[[".__S3MethodsTable__."]]
  ))
  .regFake(.ns, "fakeSame", "sameBlock", TRUE)
  .regFake(.ns, "fakeIovSame", "iovSame", TRUE)
  .regFake(.ns, "fakeFixed", "sameBlock", "fixed")
  .regFake(.ns, "fakeNone", "sameBlock", NULL)
  expect_identical(
    vapply(c("fakeSame", "fakeIovSame", "fakeFixed", "fakeNone"), .isSameBlockMethod, logical(1)),
    c(fakeSame = TRUE, fakeIovSame = TRUE, fakeFixed = TRUE, fakeNone = FALSE)
  )
  # a method that estimates the copies as one is named; one that does not
  # estimate omega ("fixed") is accepted but no alternative
  .ui <- rxode2::rxode2(.idSameMod())
  expect_identical(
    .errMsg(.nlmixr0preProcessAssertSameBlock(.ui, "saem", NULL, NULL)),
    paste0(
      "'saem' does not honour a repeated (same()) omega block\n",
      "use fo, foi, foce, focei, focep, laplace, agq, emvi, fbvi, vae, fakeIovSame or fakeSame"
    )
  )
  expect_identical(
    .errMsg(.nlmixr0preProcessAssertSameBlock(.ui, "fakeNone", NULL, NULL)),
    paste0(
      "'fakeNone' does not honour a repeated (same()) omega block\n",
      "use fo, foi, foce, focei, focep, laplace, agq, emvi, fbvi, vae, fakeIovSame or fakeSame"
    )
  )
  expect_false(.isSameBlockMethod(NA_character_))
  expect_false(.isSameBlockMethod(character(0)))
})

test_that("two occasion parameters get their occasion number, not NA", {
  skip_on_cran()
  # `.dt`'s occasion column is built from the FIRST occasion parameter's
  # `rx.<name>.<occ>` spellings; stripping the LAST one's prefix left every
  # occasion NA (and warned "NAs introduced by coercion") as soon as a
  # level carried two parameters.  Reproduces on the "theta" path too.
  .f <- .sameMod("iov.cl ~ 0.1 | occ; iov.v ~ 0.2 | occ", "cl <- exp(tcl + iov.cl)\n v <- exp(tv + iov.v)")
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.f(), .sameData(), "focei", foceiControl(print = 0, covMethod = "", iovMethod = "theta"))
  ))
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
  .f <- .sameMod("iov.cl ~ fix(0.1) | occ", "cl <- exp(tcl + iov.cl)\n v <- exp(tv)")
  .ini <- .applyIov(.f, "omega")$ui$iniDf
  .occ <- .ini[grepl("^rx\\.iov\\.cl\\.", .ini$name), ]
  expect_true(all(.occ$fix))
  expect_equal(.occ$est, rep(0.1, 2))

  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.f(), .sameData(), "focei", foceiControl(print = 0, covMethod = "", iovMethod = "omega"))
  ))
  .row <- .fit$ui$iniDf[.fit$ui$iniDf$name == "iov.cl", ]
  expect_true(.row$fix)
  expect_equal(.row$est, 0.1)
})

test_that("a fixed correlated block stays a repeated block", {
  skip_on_cran()
  # every occasion carries the same `fix` flags, so lotri still re-emits
  # `same()` -- a mismatch there would silently dissolve the repetition
  .f <- .sameMod("iov.cl + iov.v ~ fix(0.1,\n 0.03, 0.2) | occ", "cl <- exp(tcl + iov.cl)\n v <- exp(tv + iov.v)")
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
  .obj <- vapply(
    c("theta", "omega"),
    function(.m) {
      suppressWarnings(suppressMessages(
        nlmixr2est::nlmixr2(
          .f(),
          .sameData(),
          "focei",
          foceiControl(
            print = 0,
            sigdig = 7,
            covMethod = "",
            maxOuterIterations = 0L,
            maxInnerIterations = 100000L,
            iovMethod = .m
          )
        )
      ))$objf
    },
    double(1)
  )
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
    nlmixr2est::nlmixr2(.corMod(), .d, "focei", foceiControl(print = 0, covMethod = "analytic"))
  ))
  expect_false(identical(.covBaseName(.fit$covMethod), "analytic"))
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
    .uiApplyIov(rxode2::rxode2(.corMod()), "saem", .sameData(), saemControl(print = 0)),
    "not supported by 'saem'"
  )
  # asking for it outright is refused by name, not silently downgraded
  expect_error(
    .uiApplyIov(rxode2::rxode2(.corMod()), "saem", .sameData(), foceiControl(iovMethod = "omega")),
    "does not honour"
  )
})

test_that("the FOCEi family really shares the repeated block", {
  skip_on_cran()
  # Running is not enough: the reported `$occ` is read off occasion one,
  # so two independently estimated blocks would look fine.  The
  # parameter count is the tell -- 4 thetas + eta.ka + ONE 2x2 (3) = 8,
  # against 11 if each occasion were estimated on its own.
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.corMod(), .sameData(), "focei", foceiControl(print = 0, covMethod = ""))
  ))
  expect_equal(attr(logLik(.fit), "df"), 8L)
})

test_that("occasion parameters whose names share a prefix stay separate", {
  skip_on_cran()
  # The per-occasion columns are named `rx.<param>.<occ>`, and they used
  # to be selected by SUBSTRING, so `iov.v` also matched `rx.iov.v2.1`.
  # The two parameters' columns were mixed and every occasion came out
  # NA.  Reproduces on the long standing "theta" path.
  .f <- .sameMod("iov.v ~ 0.1 | occ; iov.v2 ~ 0.2 | occ", "cl <- exp(tcl + iov.v)\n v <- exp(tv + iov.v2)")
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.f(), .sameData(), "focei", foceiControl(print = 0, covMethod = "", iovMethod = "theta"))
  ))
  .occ <- .fit$iov$occ
  expect_true(all(!is.na(.occ$occ)))
  expect_equal(sort(unique(.occ$occ)), c(1L, 2L))
  # one row per subject per occasion, not doubled up
  expect_equal(nrow(.occ), 2L * length(unique(.sameData()$ID)))
  # and the two parameters are not the same column reused
  expect_false(isTRUE(all.equal(.occ$iov.v, .occ$iov.v2)))
  expect_false(any(grepl("NAs introduced", .fit$runInfo)))
})

test_that("analytic covariance bows out for a FIXED single-occasion block", {
  skip_on_cran()
  # The narrow residual of the gate: an all-`fix()`ed occasion block
  # read as "theta" when the test was on `fix` alone, and with one
  # occasion there is no repeated block for `omegaSameMap` to report.
  # "theta" is the mode that leaves these etas fixed AT ONE, so the
  # estimate has to be part of the test.
  .d <- .sameData()
  .d$occ <- 1L
  .f <- .sameMod("iov.cl + iov.v ~ fix(0.1,\n 0.03, 0.2) | occ", "cl <- exp(tcl + iov.cl)\n v <- exp(tv + iov.v)")
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.f(), .d, "focei", foceiControl(print = 0, covMethod = "analytic"))
  ))
  expect_false(identical(.covBaseName(.fit$covMethod), "analytic"))
  # the block is reported at the values it was fixed to
  expect_equal(unname(diag(.fit$omega$occ)), c(0.1, 0.2))
})

test_that("the expansion is chosen per level, not per model", {
  skip_on_cran()
  # A correlation on `occ` is no reason to change how an unrelated
  # diagonal `occ2` is expanded -- "theta" is the better conditioned
  # inner problem, so an unrelated level should keep it.
  .d <- .sameData()
  .d$occ2 <- 1 + (.d$TIME >= 9)
  .f <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- c(0, 0.7)
          eta.ka ~ 0.6
          iov.cl + iov.v ~ c(0.1,
                             0.03, 0.2) | occ
          iov.ka ~ 0.05 | occ2 })
    model({ ka <- exp(tka + eta.ka + iov.ka)
            cl <- exp(tcl + iov.cl)
            v <- exp(tv + iov.v)
            linCmt() ~ add(add.sd) })
  }
  .ini <- .uiApplyIov(rxode2::rxode2(.f()), "focei", .d, foceiControl())$ui$iniDf
  .th <- .ini[!is.na(.ini$ntheta), ]
  # the correlated level took the shared block: magnitude fixed at one
  expect_true(all(.th$fix[.th$name %in% c("iov.cl", "iov.v")]))
  expect_equal(.th$est[.th$name == "iov.cl"], 1)
  # the diagonal level kept the old expansion: a free SD-scale magnitude
  expect_false(.th$fix[.th$name == "iov.ka"])
  expect_equal(.th$est[.th$name == "iov.ka"], sqrt(0.05))
  # ... and its per-occasion etas are still unit variance and fixed
  .ka <- .ini[grepl("^rx\\.iov\\.ka\\.", .ini$name), ]
  expect_true(all(.ka$fix))
  expect_equal(.ka$est, rep(1, 2))
  expect_true(all(.ka$condition == "id"))
  # only the correlated level is repeated
  expect_equal(sum(grepl(":same:", .ini$condition, fixed = TRUE)), 3L)
})

test_that("a mixed correlated/diagonal model fits and reports both levels", {
  skip_on_cran()
  .d <- .sameData()
  .d$occ2 <- 1 + (.d$TIME >= 9)
  .f <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- c(0, 0.7)
          eta.ka ~ 0.6
          iov.cl + iov.v ~ c(0.1,
                             0.03, 0.2) | occ
          iov.ka ~ 0.05 | occ2 })
    model({ ka <- exp(tka + eta.ka + iov.ka)
            cl <- exp(tcl + iov.cl)
            v <- exp(tv + iov.v)
            linCmt() ~ add(add.sd) })
  }
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.f(), .d, "focei", foceiControl(print = 0, covMethod = ""))
  ))
  expect_equal(dim(.fit$omega$occ), c(2L, 2L))
  expect_equal(dim(.fit$omega$occ2), c(1L, 1L))
  .ini <- .fit$ui$iniDf
  expect_equal(
    sort(.ini$name[.ini$condition %in% "occ" & !is.na(.ini$condition)]),
    sort(c("iov.cl", "iov.v", "(iov.cl,iov.v)"))
  )
  expect_equal(.ini$name[.ini$condition %in% "occ2" & !is.na(.ini$condition)], "iov.ka")
  expect_true(all(!is.na(.fit$iov$occ$occ)))
  expect_true(all(!is.na(.fit$iov$occ2$occ2)))
})

test_that("foceiFixed() stays aligned when a fixed eta follows a copy", {
  skip_on_cran()
  # `rxUiGet.foceiFixed()` returns `c(theta flags, eta-row flags)` and is
  # consumed POSITIONALLY against the optimizer's `[theta | omegaTheta]`
  # layout.  A `same()` copy contributes no omega theta of its own, so
  # its rows must not appear -- otherwise every flag after them is off by
  # one and the WRONG parameter is held fixed, in a fit that still
  # converges.  It only bites when a fixed eta row sits AFTER a copy row.
  .d <- .sameData()
  .d$occ2 <- 1 + (.d$TIME >= 9)
  .f <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- c(0, 0.7)
          eta.ka ~ 0.6
          iov.cl + iov.v ~ c(0.1,
                             0.03, 0.2) | occ
          iov.ka ~ fix(0.05) | occ2 })
    model({ ka <- exp(tka + eta.ka + iov.ka)
            cl <- exp(tcl + iov.cl)
            v <- exp(tv + iov.v)
            linCmt() ~ add(add.sd) })
  }
  .ui <- .uiApplyIov(rxode2::rxode2(.f()), "focei", .d, foceiControl())$ui
  .fx <- rxUiGet.foceiFixed(list(.ui))
  .nth <- sum(!is.na(.ui$iniDf$ntheta))
  .rxInv <- rxode2::rxSymInvCholCreate(.ui$omega, "sqrt", same = .ui$omegaSameMap)
  # the contract: one flag per theta, then one per OMEGA PARAMETER
  expect_equal(length(.fx), .nth + length(.rxInv$theta))
  # three magnitude thetas are fixed: iov.cl and iov.v at one (shared
  # block) and iov.ka by the user's own fix()
  expect_equal(sum(.fx[seq_len(.nth)]), 3L)
  # The omega half is eta.ka, then the ONE 2x2 block (3), then iov.ka's
  # two occasion etas -- which "theta" pins at unit variance.  The
  # alignment is the point: the two TRUEs must be the LAST two entries.
  # Unfiltered, the copy rows would push them two places along and hold
  # the wrong parameters fixed in a fit that still converges.
  .om <- .fx[-seq_len(.nth)]
  expect_equal(length(.om), 6L)
  expect_equal(which(.om), c(5L, 6L))
})

test_that("three parameters over three occasions repeat correctly", {
  skip_on_cran()
  .d <- .sameData()
  .d$occ <- 1 + (.d$TIME >= 3) + (.d$TIME >= 7)
  .f <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- c(0, 0.7)
          eta.ka ~ 0.6
          iov.ka + iov.cl + iov.v ~ c(0.1,
                                      0.02, 0.2,
                                      0.01, 0.03, 0.3) | occ })
    model({ ka <- exp(tka + eta.ka + iov.ka)
            cl <- exp(tcl + iov.cl)
            v <- exp(tv + iov.v)
            linCmt() ~ add(add.sd) })
  }
  .ui <- .uiApplyIov(rxode2::rxode2(.f()), "focei", .d, foceiControl(iovMethod = "omega"))$ui
  .om <- .ui$omega
  expect_equal(dim(.om), c(10L, 10L))
  # eta.ka, then three identical 3x3 blocks
  expect_equal(unname(.om[2:4, 2:4]), unname(.om[5:7, 5:7]))
  expect_equal(unname(.om[2:4, 2:4]), unname(.om[8:10, 8:10]))
  # each occasion mirrors occasion ONE, not the occasion before it
  expect_equal(.ui$omegaSameMap, c(0L, 0L, 0L, 0L, 2L, 3L, 4L, 2L, 3L, 4L))
  # 1 (eta.ka) + 6 (one 3x3 block) rather than 1 + 18
  expect_equal(
    length(
      rxode2::rxSymInvCholCreate(
        .om,
        "sqrt",
        same = .ui$omegaSameMap
      )$theta
    ),
    7L
  )
})

test_that("a repeated occasion-level block stops before the fit, naming its etas", {
  skip_on_cran()
  .msg <- "same() is not supported for an occasion-level block: "
  .one <- .sameMod(
    "iov.cl ~ 0.1 | occ\n iov.v ~ same() | occ",
    "cl <- exp(tcl + iov.cl)\n v <- exp(tv + iov.v)"
  )
  # every method that handles IOV, saem's own handling included, whether or
  # not it honours a repeated subject-level block
  for (.est in c("focei", "saem", "imp", "laplace", "agq", "npag", "vae")) {
    .e <- expect_error(nlmixr2(.one, .sameData(), .est, control = list(print = 0)), label = .est)
    expect_identical(conditionMessage(.e), paste0(.msg, "'iov.v' | occ"), label = .est)
  }
  .e <- expect_error(nlmixr2(.one, .sameData(), "saem", control = saemControl(print = 0, iovMethod = "theta")))
  expect_identical(conditionMessage(.e), paste0(.msg, "'iov.v' | occ"))
  # a repeated 2x2 block names both etas
  .two <- .sameMod(
    "iov.cl + iov.v ~ c(0.1, 0.01, 0.1) | occ\n iov.ka + iov.q ~ same() | occ",
    "ka <- exp(tka + eta.ka + iov.ka)\n cl <- exp(tcl + iov.cl)\n v <- exp(tv + iov.v + iov.q)"
  )
  .e <- expect_error(nlmixr2(.two, .sameData(), "focei", control = list(print = 0)))
  expect_identical(conditionMessage(.e), paste0(.msg, "'iov.ka', 'iov.q' | occ"))
  # one group per occasion column
  .lvl <- .sameMod(
    "iov.cl ~ 0.1 | occ\n iov.v ~ same() | occ\n iov.ka ~ 0.1 | visit\n iov.q ~ same() | visit",
    "ka <- exp(tka + eta.ka + iov.ka)\n cl <- exp(tcl + iov.cl)\n v <- exp(tv + iov.v + iov.q)"
  )
  .e <- expect_error(.iovSameRefuse(rxode2::rxode2(.lvl()), "focei"))
  expect_identical(conditionMessage(.e), paste0(.msg, "'iov.v' | occ; 'iov.q' | visit"))
})

test_that("the same() refusal leaves id-level blocks and simulation alone", {
  skip_on_cran()
  # eta.v repeats eta.cl at the subject level, beside an occasion eta
  .id <- .sameMod(
    "eta.cl ~ 0.3\n eta.v ~ same()\n iov.cl ~ 0.1 | occ",
    "cl <- exp(tcl + eta.cl + iov.cl)\n v <- exp(tv + eta.v)"
  )
  expect_null(.iovSameRefuse(rxode2::rxode2(.id()), "focei"))
  # rxSolve simulates the occasion block itself
  .one <- .sameMod(
    "iov.cl ~ 0.1 | occ\n iov.v ~ same() | occ",
    "cl <- exp(tcl + iov.cl)\n v <- exp(tv + iov.v)"
  )
  expect_null(.iovSameRefuse(rxode2::rxode2(.one()), "rxSolve"))
  .s <- suppressMessages(nlmixr2(.one, .sameData(), "rxSolve"))
  expect_equal(nrow(as.data.frame(.s)), sum(.sameData()$EVID == 0))
})
