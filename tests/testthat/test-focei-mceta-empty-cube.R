nmTest({
  # Regression test for a Cube::slice() out-of-bounds abort in the FOCEi inner
  # loop (src/inner.cpp `innerOpt1`).  When `mceta >= 1` the inner loop reads
  # pre-drawn ETA samples from `op_focei.mcetaSamples.slice(id)`, but that cube
  # is only filled when `maxInnerIterations > 0`.  A pure evaluation -- e.g. the
  # covariance step, or `nlmixr2extra::linearize()` -- runs with
  # `maxInnerIterations == 0`, leaving the cube empty (n_slices == 0); the read
  # then indexed an empty cube and threw `std::out_of_range`
  # ("Cube::slice(): index out of bounds") from inside the OpenMP parallel
  # region, which is uncaught and aborts R.  inner.cpp now guards the read with
  # `id < n_slices`.  This fit reproduces the condition directly (mceta >= 1 with
  # maxInnerIterations == 0); before the fix it aborted the R process.
  test_that("FOCEi mceta does not index an empty mcetaSamples cube (maxInnerIterations=0)", {
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    fit <- suppressWarnings(suppressMessages(
      nlmixr(one.cmt, nlmixr2data::theo_sd, "focei",
             foceiControl(mceta = 10, maxInnerIterations = 0))))
    expect_true(inherits(fit, "nlmixr2FitData"))
  })
})
