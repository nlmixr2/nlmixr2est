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

  test_that("etaMat drops the mixture bookkeeping columns", {
    ## $eta carries ID, and for a mixture fit nmObjGet.ranef merges in a mixnum
    ## column.  Neither is an eta.  Letting either reach an etaMat trips
    ## foceiSetup_'s column check ("The etaMat must have the same number of ETAs
    ## (cols) as the model"), which is what broke every mixture refit -- $cov,
    ## addCwres, .setOfvFo and nlmixr2(fit, ...) all round-trip fit$etaMat.
    .eta <- data.frame(ID = 1:3, eta.ka = c(0.1, -0.2, 0.3),
                       eta.cl = c(-0.1, 0.2, 0.0), mixnum = c(1L, 2L, 1L))
    expect_equal(colnames(.nmDropNonEtaCols(.eta)), c("eta.ka", "eta.cl"))
    ## the accessor itself, with and without the mixture column
    expect_equal(colnames(nmObjGet.etaMat(list(list(eta = .eta, iov = NULL)))),
                 c("eta.ka", "eta.cl"))
    expect_equal(colnames(nmObjGet.etaMat(list(list(eta = .eta[, 1:3], iov = NULL)))),
                 c("eta.ka", "eta.cl"))
    ## MIXEST is the same bookkeeping under its focei spelling
    .eta2 <- .eta
    names(.eta2)[4] <- "MIXEST"
    expect_equal(colnames(.nmDropNonEtaCols(.eta2)), c("eta.ka", "eta.cl"))
    ## an ordinary eta whose name merely contains "mix" is NOT dropped
    .eta3 <- data.frame(ID = 1:2, eta.mixup = c(0.1, 0.2))
    expect_equal(colnames(.nmDropNonEtaCols(.eta3)), "eta.mixup")
  })
})
