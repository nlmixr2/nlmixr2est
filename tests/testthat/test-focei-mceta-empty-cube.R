nmTest({
  # Regression test: mceta >= 1 with maxInnerIterations == 0 left
  # mcetaSamples empty, so innerOpt1() indexed an out-of-bounds cube slice
  # and aborted R (src/inner.cpp).
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
