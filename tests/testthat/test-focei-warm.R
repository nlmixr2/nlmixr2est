nmTest({
  test_that("foceiControl(warm=) option mapping", {
    expect_equal(foceiControl()$warm, 1L)
    expect_equal(foceiControl(warm="calc")$warm, 1L)
    expect_equal(foceiControl(warm="save")$warm, 0L)
    expect_equal(foceiControl(warm=0L)$warm, 0L)
    expect_equal(foceiControl(warm=1L)$warm, 1L)
    expect_error(foceiControl(warm="bogus"))
    .ctl <- foceiControl(warm="save")
    expect_equal(do.call(foceiControl, .ctl)$warm, 0L)
  })

  test_that("focei warm='calc' matches warm='save'", {
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- log(c(0, 2.7, 100))
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

    .fit <- function(warm, maxOuterIterations) {
      suppressWarnings(suppressMessages(
        nlmixr(one.cmt, nlmixr2data::theo_sd, "focei",
               foceiControl(maxOuterIterations=maxOuterIterations,
                            covMethod="", calcTables=FALSE, print=0,
                            warm=warm))))
    }

    # posthoc: same inner problems converged to the same etas/objective
    .p1 <- .fit("calc", 0L)
    .p2 <- .fit("save", 0L)
    expect_equal(.p1$objf, .p2$objf, tolerance=1e-4)
    expect_equal(as.data.frame(.p1$eta), as.data.frame(.p2$eta),
                 tolerance=1e-4)

    # short optimization run finishes and agrees
    .f1 <- .fit("calc", 5L)
    .f2 <- .fit("save", 5L)
    expect_true(inherits(.f1, "nlmixr2FitCore"))
    expect_true(inherits(.f2, "nlmixr2FitCore"))
    expect_equal(.f1$objf, .f2$objf, tolerance=1e-2)
  })
})
