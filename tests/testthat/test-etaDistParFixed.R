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
    # the dist() parameters carry a back-transformation now.  They are not
    # mu-referenced -- they only ever reach the model inside the declaration --
    # so nothing recorded their scale and parFixed printed them raw: a
    # clearance of 1.5 showed as 0.405.
    for (.p in c("lclm", "lclrv", "lv1m", "lv1rv")) {
      expect_equal(unname(.f$parFixedDf[.p, "Back-transformed"]),
                   exp(unname(.f$parFixedDf[.p, "Estimate"])), info = .p)
    }
    # at maxOuterIterations = 0 those are the values the user typed
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
})
