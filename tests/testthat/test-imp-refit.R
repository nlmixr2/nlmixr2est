# Re-fitting a completed imp-family fit from the fit object,
# nlmixr2(fit, est = ...).  This was broken outright: getValidNlmixrCtl.impmap
# round-trips through impmapControl(), which forwards anything it does not
# recognise to foceiControl(), and the four per-model M-step index maps
# .impmapFamilyFit stamps on the runtime control are arguments of neither.  So
# every re-fit died with "unused argument: 'impMuThetaIdx', ...".
#
# The second half is the rule that only becomes reachable once the first is
# fixed: est="imp" stamps mapIter = 0 on its control (that is what imp means),
# so adopting an imp fit's control under est="impmap" would silently run a
# method that never re-optimizes the mode.  `est` has to win -- but only over
# an INHERITED 0, never over one the user asked for.
nmTest({
  .one <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  .ctl <- function(...) {
    impmapControl(print = 0L, nIter = 3L, isample = 200L, covMethod = "",
                  calcTables = FALSE, ...)
  }

  test_that("a completed fit's own control survives re-validation", {
    # the pure control-surface half, no re-fit: do.call(impmapControl, ctl) is
    # exactly what getValidNlmixrCtl.impmap does
    .c0 <- .ctl()
    # simulate what .impmapFamilyFit stamps on the runtime control
    .c0$impMuThetaIdx <- 0:1
    .c0$impMuEtaIdx <- 0:1
    .c0$impThetaSensIdx <- integer(0)
    .c0$impOmegaFixedEta <- integer(0)
    .c1 <- do.call(impmapControl, .c0)
    expect_s3_class(.c1, "impmapControl")
    # carried through, so the round-trip is idempotent
    for (.nm in .impmapIdxMapNames) {
      expect_identical(.c1[[.nm]], .c0[[.nm]])
    }
    expect_identical(do.call(impmapControl, .c1)[[.impmapIdxMapNames[1]]],
                     .c0[[.impmapIdxMapNames[1]]])
    # and they are still stripped when down-converting to a plain foceiControl
    expect_true(all(.impmapIdxMapNames %in% .impmapIsControlNames))
  })

  test_that("every imp-family fit can be re-fit from the fit object", {
    .dat <- nlmixr2data::theo_sd
    .f <- suppressWarnings(suppressMessages(
      nlmixr2(.one, .dat, "imp", impControl(print = 0L, nIter = 3L,
                                            isample = 200L, covMethod = "",
                                            calcTables = FALSE))))
    for (.e in c("imp", "impmap", "qrpem")) {
      .r <- suppressWarnings(suppressMessages(nlmixr2(.f, est = .e)))
      expect_true(is.finite(.r$objf))
      expect_true(all(.impmapIdxMapNames %in% names(.r$env$control)))
    }
  })

  test_that("est wins over a mapIter inherited from an imp fit", {
    .dat <- nlmixr2data::theo_sd
    .fi <- suppressWarnings(suppressMessages(
      nlmixr2(.one, .dat, "imp", impControl(print = 0L, nIter = 3L,
                                            isample = 200L, covMethod = "",
                                            calcTables = FALSE))))
    expect_identical(.fi$env$impMapIter, 0L)
    # re-fitting as impmap must NOT inherit "never re-center"
    .r <- suppressWarnings(suppressMessages(nlmixr2(.fi, est = "impmap")))
    expect_identical(.r$env$impMapIter, 1L)
    .rq <- suppressWarnings(suppressMessages(nlmixr2(.fi, est = "qrpem")))
    expect_identical(.rq$env$impMapIter, 1L)
    # ... and re-fitting as imp still means imp
    .ri <- suppressWarnings(suppressMessages(nlmixr2(.fi, est = "imp")))
    expect_identical(.ri$env$impMapIter, 0L)
  })

  test_that("est does not override a mapIter the user asked for", {
    .dat <- nlmixr2data::theo_sd
    # a user's own impmapControl(mapIter = 0) carries no `est`, which is what
    # distinguishes it from an inherited one
    .c <- .ctl(mapIter = 0L)
    expect_null(.c$est)
    expect_identical(do.call(impmapControl, .c)$mapIter, 0L)

    .fu <- suppressWarnings(suppressMessages(nlmixr2(.one, .dat, "impmap", .c)))
    expect_identical(.fu$env$impMapIter, 0L)
    # re-fitting an IMPMAP fit as impmap preserves it -- only an est="imp"
    # source triggers the reset
    .r <- suppressWarnings(suppressMessages(nlmixr2(.fu, est = "impmap")))
    expect_identical(.r$env$impMapIter, 0L)
  })
})
