nmTest({
  test_that("etaMat works", {

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
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

    # This test only inspects how etaMat is stored in foceiControl, not any
    # numeric result, so cap the outer iterations to keep the fits fast.
    f <- .nlmixr(one.cmt, theo_sd, "focei",
                 foceiControl(maxOuterIterations=0L))

    expect_null(f$foceiControl$etaMat)

    f2 <- .nlmixr(f, est="focei", control=foceiControl(outerOpt="bobyqa", maxOuterIterations=0L))
    expect_true(inherits(f2$foceiControl$etaMat, "matrix"))

    f2 <- .nlmixr(f, "focei", foceiControl(outerOpt="bobyqa", maxOuterIterations=0L))
    expect_true(inherits(f2$foceiControl$etaMat, "matrix"))

    f3 <- .nlmixr(f, est="focei", control=foceiControl(outerOpt="bobyqa", etaMat=f, maxOuterIterations=0L))
    expect_true(inherits(f3$foceiControl$etaMat, "matrix"))

    f4 <- .nlmixr(f, est="focei", control=foceiControl(outerOpt="bobyqa",
                                                       etaMat=NA, maxOuterIterations=0L))
    expect_true(is.na(f4$foceiControl$etaMat))

    f4 <- .nlmixr(f, "focei", foceiControl(outerOpt="bobyqa",
                                           etaMat=NA, maxOuterIterations=0L))

    expect_true(is.na(f4$foceiControl$etaMat))

    f4 <- .nlmixr(f, foceiControl(outerOpt="bobyqa",
                                  etaMat=NA, maxOuterIterations=0L))

    expect_true(is.na(f4$foceiControl$etaMat))

  })
})
