# Re-fitting a completed imp-family fit from the fit object,
# nlmixr2(fit, est = ...).  This was broken outright: getValidNlmixrCtl.impmap
# round-trips through impmapControl(), which forwards anything it does not
# recognise to foceiControl(), and the four per-model M-step index maps
# .impmapFamilyFit stamps on the runtime control are arguments of neither.  So
# every re-fit died with "unused argument: 'impMuThetaIdx', ...".
#
# The second half is the rule that only becomes reachable once the first is
# fixed.  Several fields on a completed fit's control were put there by that
# fit's est rather than by the user: est="imp" stamps mapIter = 0 (that IS what
# imp means) and est="qrpem" stamps qr = TRUE, sir = TRUE.  Carrying either into
# a different method silently runs a different algorithm than the one asked for
# -- a re-fit as qrpem inheriting an imp control drew plain Monte-Carlo samples
# and was still labelled QRPEM.  `est` has to win, but ONLY over a value another
# est stamped, never over one the user wrote themselves (which is why the rule
# is keyed on the control's `est` field, not on the value).
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

  test_that("est wins over qr/sir inherited from a qrpem fit, both ways", {
    .dat <- nlmixr2data::theo_sd
    .fi <- suppressWarnings(suppressMessages(
      nlmixr2(.one, .dat, "imp", impControl(print = 0L, nIter = 2L,
                                            isample = 150L, covMethod = "",
                                            calcTables = FALSE))))
    expect_false(.fi$env$impQr)
    # est="qrpem" IS impmapControl(qr=TRUE, sir=TRUE).  Re-fitting an imp fit
    # as qrpem must not draw plain Monte-Carlo samples and call it QRPEM.
    .rq <- suppressWarnings(suppressMessages(nlmixr2(.fi, est = "qrpem")))
    expect_true(.rq$env$impQr)
    expect_true(.rq$env$impSir)

    # and the other direction: qrpem's qr/sir must not leak into imp/impmap
    .fq <- suppressWarnings(suppressMessages(
      nlmixr2(.one, .dat, "qrpem", qrpemControl(print = 0L, nIter = 2L,
                                                isample = 150L, covMethod = "",
                                                calcTables = FALSE))))
    expect_true(.fq$env$impQr)
    .ri <- suppressWarnings(suppressMessages(nlmixr2(.fq, est = "imp")))
    expect_false(.ri$env$impQr)
    expect_false(.ri$env$impSir)
  })

  test_that("est does not override a qr/sir the user asked for", {
    # a user-built control carries no `est`, which is what distinguishes it
    # from one inherited off a completed fit
    .c <- qrpemControl(qr = FALSE, sir = FALSE)
    expect_null(.c$est)
    .v <- getValidNlmixrCtl.qrpem(structure(list(.c), class = "qrpem"))
    expect_false(.v$qr)
    expect_false(.v$sir)
    # ... and an impmapControl(qr=TRUE) survives validation as impmap
    .c2 <- impmapControl(qr = TRUE)
    .v2 <- getValidNlmixrCtl.impmap(structure(list(.c2), class = "impmap"))
    expect_true(.v2$qr)
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
