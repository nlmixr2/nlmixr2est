# impmapControl(nBurn=, burnFreezeOmega=) -- burn-in EM iterations that let the
# gamma / auto controllers settle before the estimates they influence are
# judged.  The mechanism assertions are (a) the iteration budget grows by nBurn
# rather than being carved out of nIter, (b) Omega literally does not move for
# the first nBurn iterations under burnFreezeOmega while the thetas do, and
# (c) nBurn = 0 leaves the fit bit-identical.  Two small fits; essential subset.
nmTest({
  test_that("nBurn / burnFreezeOmega validation and round-trip", {
    expect_identical(impmapControl()$nBurn, 0L)
    expect_false(impmapControl()$burnFreezeOmega)
    expect_identical(impmapControl(nBurn = 5L)$nBurn, 5L)
    expect_true(impmapControl(burnFreezeOmega = TRUE)$burnFreezeOmega)
    # inherited by the imp / qrpem shims
    expect_identical(impControl(nBurn = 4L)$nBurn, 4L)
    expect_identical(qrpemControl(nBurn = 2L)$nBurn, 2L)

    expect_error(impmapControl(nBurn = -1L), "nBurn")
    expect_error(impmapControl(nBurn = c(1L, 2L)), "nBurn")
    expect_error(impmapControl(nBurn = NA_integer_), "nBurn")
    expect_error(impmapControl(burnFreezeOmega = NA), "burnFreezeOmega")
    expect_error(impmapControl(burnFreezeOmega = c(TRUE, FALSE)), "burnFreezeOmega")

    .ctl <- impmapControl(nBurn = 3L, burnFreezeOmega = TRUE)
    expect_identical(do.call(impmapControl, .ctl)$nBurn, 3L)
    expect_true(do.call(impmapControl, .ctl)$burnFreezeOmega)

    # stripped when down-converting to a plain foceiControl, and declared inert
    # for the nonparametric engines
    expect_true(all(c("nBurn", "burnFreezeOmega") %in% .impmapIsControlNames))
    expect_true(all(c("nBurn", "burnFreezeOmega") %in% .npInertImpCtl))
  })

  test_that("burn-in iterations are extra, and freeze Omega without freezing theta", {
    one.cmt <- function() {
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
    .dat <- nlmixr2data::theo_sd
    .run <- function(...) {
      suppressWarnings(suppressMessages(
        nlmixr2(one.cmt, .dat, "impmap",
                impmapControl(print = 0L, nIter = 5L, isample = 100L,
                              covMethod = "", calcTables = FALSE, ...))))
    }
    .f0 <- .run()
    .fb <- .run(nBurn = 3L, burnFreezeOmega = TRUE)

    # reported back so a reader can split the traces
    expect_identical(.f0$env$impNburn, 0L)
    expect_identical(.fb$env$impNburn, 3L)
    expect_false(.f0$env$impBurnFreezeOmega)
    expect_true(.fb$env$impBurnFreezeOmega)

    # EXTRA, not carved out of nIter: nBurn=3 + nIter=5 runs 8 iterations
    expect_identical(.f0$env$impIter, 5L)
    expect_identical(.fb$env$impIter, 8L)
    expect_length(.fb$env$impObjTrace, 8L)

    # the frozen half: the omega entries of $parHist do not move across the
    # first nBurn iterations, and then do afterwards
    .ph <- .fb$env$parHistData
    .ph <- .ph[.ph$type == "Scaled", , drop = FALSE]
    .om <- grep("^o[0-9]+$", names(.ph), value = TRUE)
    expect_true(length(.om) > 0L)
    .omBurn <- .ph[seq_len(3L), .om, drop = FALSE]
    for (.j in .om) expect_equal(length(unique(.omBurn[[.j]])), 1L)
    # ... and unfreezing actually releases it
    expect_false(isTRUE(all.equal(unlist(.ph[4L, .om]), unlist(.ph[3L, .om]),
                                  tolerance = 1e-8)))

    # the thetas are NOT frozen -- the whole M-step still ran
    expect_false(isTRUE(all.equal(.ph$tka[3L], .ph$tka[1L], tolerance = 1e-8)))
  })

  test_that("convergence cannot fire until the window clears the burn-in", {
    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .dat <- nlmixr2data::theo_sd
    # a window of 2 with a loose ctol converges almost immediately without a
    # burn-in; with one, it cannot stop before iteration nBurn + nConvWindow
    .fb <- suppressWarnings(suppressMessages(
      nlmixr2(one.cmt, .dat, "impmap",
              impmapControl(print = 0L, nIter = 20L, isample = 100L,
                            nBurn = 4L, nConvWindow = 2L, ctol = 1e-2,
                            gammaRule = "floor", covMethod = "",
                            calcTables = FALSE))))
    expect_gte(.fb$env$impIter, 4L + 2L + 1L)
  })
})
