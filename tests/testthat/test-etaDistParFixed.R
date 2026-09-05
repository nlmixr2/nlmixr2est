nmTest({
  test_that("parFixed reports the expansion the way the model was written", {
    # rxEtaDistExpand() leaves two kinds of row a user never wrote:
    #
    #  * the latent standard normals (rxz.<eta>), whose variance is fixed at
    #    one by construction -- they are not estimates and print as NA.
    #  * the copula correlations (rxCor.<i>.<j>), which ARE estimated, but the
    #    number carried is the UNCONSTRAINED parameter -- the expansion writes
    #    tanh() around it so the optimizer can range over the real line.
    #
    # Printed raw next to real estimates the second reads as a correlation and
    # overstates it: 0.549 is a correlation of 0.5, and a fitted 1.047 is 0.78,
    # which is not even in the legal range for one.  Two summaries written
    # while developing this feature quoted that raw number as a correlation.
    .m <- function() {
      ini({
        lclm <- log(1.5)
        lclrv <- log(0.1)
        lv1m <- log(4.7)
        lv1rv <- log(0.09)
        tka <- 0.45
        eta.cl + eta.v1 ~ c(1,
                            0.5, 1)
        dist(eta.cl) ~ dgamma(shape = 1 / exp(lclrv),
                              rate = 1 / (exp(lclrv) * exp(lclm)))
        dist(eta.v1) ~ dgamma(shape = 1 / exp(lv1rv),
                              rate = 1 / (exp(lv1rv) * exp(lv1m)))
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- eta.cl
        v <- eta.v1
        linCmt() ~ add(add.sd)
      })
    }
    .f <- suppressWarnings(nlmixr2(.m(), nlmixr2data::theo_sd, est = "focei",
                                   control = foceiControl(print = 0,
                                                          covMethod = "",
                                                          maxOuterIterations = 0)))
    # the latents are gone from both the frame and the printed copy
    expect_false(any(grepl("^rxz[.]", rownames(.f$parFixedDf))))
    expect_false(any(grepl("^rxz[.]", rownames(.f$parFixed))))
    # ...and the parameters the user DID write are still all there
    for (.p in c("lclm", "lclrv", "lv1m", "lv1rv", "tka", "add.sd")) {
      expect_true(.p %in% rownames(.f$parFixedDf), info = .p)
    }
    # The dist() parameters and the copula correlation are BOTH already
    # reported on their natural scale, and were before any of this work:
    # rxode2 records curEval = "exp" for a theta appearing inside exp() even
    # when it is not mu-referenced, and the expansion sets
    # backTransform = "tanh" on the correlation rows.  These guard against
    # either being lost -- they are not describing something that was fixed.
    # (Both were "fixed" here at one point, on the strength of my reading the
    # raw Estimate column instead of the Back-transformed one.)
    for (.p in c("lclm", "lclrv", "lv1m", "lv1rv")) {
      expect_equal(unname(.f$parFixedDf[.p, "Back-transformed"]),
                   exp(unname(.f$parFixedDf[.p, "Estimate"])), info = .p)
    }
    expect_equal(unname(.f$parFixedDf["lclm", "Back-transformed"]), 1.5)
    expect_equal(unname(.f$parFixedDf["lclrv", "Back-transformed"]), 0.1)

    # the copula correlation was ALREADY reported correctly -- the expansion
    # sets backTransform = "tanh" -- so this only guards against that being
    # lost.  (It was briefly "fixed" here on the strength of a misreading of
    # the raw Estimate column.)
    .c <- grep("^rxCor[.]", rownames(.f$parFixedDf))
    expect_equal(length(.c), 1L)
    .raw <- unname(.f$parFixedDf[.c, "Estimate"])
    .bck <- unname(.f$parFixedDf[.c, "Back-transformed"])
    expect_equal(.bck, tanh(.raw))
    expect_equal(.bck, 0.5)
    expect_true(abs(.bck) <= 1)
  })

  test_that("the declared block is reported in $omega, and etaDistCor works", {
    # The user writes `eta.cl + eta.v1 ~ c(1, 0.5, 1)` -- a covariance block
    # with the correlation in it.  What $omega used to return was the
    # expansion's internals: a 2x2 identity named rxz.eta.cl / rxz.eta.v1, with
    # the fitted correlation off in a rxCor.* theta.  The block the user wrote
    # was nowhere in the output.
    #
    # The latents are standard normals, so their covariance matrix IS the
    # correlation matrix (unit diagonal is what the declaration requires), and
    # the block goes back into $omega on the covariance scale under the names
    # the model used.  $omegaR then derives the correlation view through the
    # machinery every other model uses.
    #
    # etaDistCor returned NULL for this same model: `etaDistInfo` is recorded
    # on the ui rxEtaDistExpand() returns and does not survive onto fit$ui, so
    # both now recover the blocks from the fit itself.
    .m <- function() {
      ini({
        lclm <- log(1.5)
        lclrv <- log(0.1)
        lv1m <- log(4.7)
        lv1rv <- log(0.09)
        tka <- 0.45
        eta.cl + eta.v1 ~ c(1,
                            0.5, 1)
        dist(eta.cl) ~ dgamma(shape = 1 / exp(lclrv),
                              rate = 1 / (exp(lclrv) * exp(lclm)))
        dist(eta.v1) ~ dgamma(shape = 1 / exp(lv1rv),
                              rate = 1 / (exp(lv1rv) * exp(lv1m)))
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- eta.cl
        v <- eta.v1
        linCmt() ~ add(add.sd)
      })
    }
    .f <- suppressWarnings(nlmixr2(.m(), nlmixr2data::theo_sd, est = "focei",
                                   control = foceiControl(print = 0,
                                                          covMethod = "",
                                                          maxOuterIterations = 0)))
    # named as the model named them, not as the expansion did
    expect_setequal(colnames(.f$omega), c("eta.cl", "eta.v1"))
    expect_false(any(grepl("^rxz[.]", colnames(.f$omega))))
    # unit diagonal, and the declared correlation off it
    expect_equal(unname(diag(.f$omega)), c(1, 1))
    expect_equal(unname(.f$omega["eta.cl", "eta.v1"]), 0.5)
    expect_equal(unname(.f$omega["eta.v1", "eta.cl"]), 0.5)
    # so the standard correlation accessor works on it
    expect_equal(unname(.f$omegaR["eta.cl", "eta.v1"]), 0.5)
    # and etaDistCor returns the block rather than NULL
    .c <- .f$etaDistCor
    expect_false(is.null(.c))
    expect_equal(length(.c), 1L)
    expect_equal(unname(.c[[1]]["eta.cl", "eta.v1"]), 0.5)
    # the uncertainty stays in the output: moving the value into $omega must
    # not take the standard error out with it
    expect_true("rxCor.eta.v1.eta.cl" %in% rownames(.f$parFixedDf))
  })
})
