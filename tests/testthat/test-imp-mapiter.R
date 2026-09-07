# impmapControl(mapIter=) -- the MAP-assist period.  Before this it was accepted
# by the control and never read by the kernel (grep mapIter src/ found nothing),
# so a user raising it changed nothing.  These tests assert the MECHANISM: the
# objective trace must diverge at exactly the iteration the period skips, not
# merely that the fits look similar.  Two tiny fits, kept in the essential
# (non-slow) subset.
nmTest({
  test_that("mapIter validation and round-trip", {
    expect_identical(impmapControl()$mapIter, 1L)
    expect_identical(impmapControl(mapIter = 4L)$mapIter, 4L)
    expect_identical(impmapControl(mapIter = 0L)$mapIter, 0L)
    # est="imp" is the no-MAP variant and pins the period to 0
    expect_identical(impControl()$mapIter, 0L)
    expect_identical(impControl(mapIter = 5L)$mapIter, 0L)
    # qrpem inherits the impmap default
    expect_identical(qrpemControl()$mapIter, 1L)

    expect_error(impmapControl(mapIter = -1L), "mapIter")
    expect_error(impmapControl(mapIter = c(1L, 2L)), "mapIter")
    expect_error(impmapControl(mapIter = NA_integer_), "mapIter")

    # a control round-tripped through itself is idempotent
    .ctl <- impmapControl(mapIter = 3L)
    expect_identical(do.call(impmapControl, .ctl)$mapIter, 3L)

    # the name is stripped when down-converting to a plain foceiControl
    expect_true("mapIter" %in% .impmapIsControlNames)
  })

  test_that("mapIter skips the MAP re-centering on exactly the right iterations", {
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
    .run <- function(mapIter) {
      suppressWarnings(suppressMessages(
        nlmixr2(one.cmt, .dat, "impmap",
                impmapControl(print = 0L, nIter = 5L, isample = 100L,
                              mapIter = mapIter, covMethod = "",
                              calcTables = FALSE))))
    }
    .f1 <- .run(1L)
    .f3 <- .run(3L)

    # the period is reported back, so a reader can tell what actually ran
    expect_identical(.f1$env$impMapIter, 1L)
    expect_identical(.f3$env$impMapIter, 3L)

    .t1 <- .f1$env$impObjTrace
    .t3 <- .f3$env$impObjTrace
    expect_length(.t1, 5L)
    expect_length(.t3, 5L)

    # iteration 0 never re-centers under either setting, so it is identical
    expect_equal(.t3[1], .t1[1], tolerance = 1e-12)
    # iteration 1 DOES re-center at mapIter=1 and does NOT at mapIter=3 -- this
    # is the assertion that the kernel actually reads the period
    expect_false(isTRUE(all.equal(.t3[2], .t1[2], tolerance = 1e-8)))
    # ... and mapIter=3 re-centers again at iteration 3, so it tracks a fit that
    # never re-centers only up to there
    .t0 <- .run(0L)$env$impObjTrace
    expect_equal(.t3[1:3], .t0[1:3], tolerance = 1e-12)
    expect_false(isTRUE(all.equal(.t3[4], .t0[4], tolerance = 1e-8)))

    # a skipped MAP is a cost/accuracy trade, not a different problem: both
    # still land on the same objective scale
    expect_equal(.f3$objf, .f1$objf, tolerance = 0.05)
  })
})
