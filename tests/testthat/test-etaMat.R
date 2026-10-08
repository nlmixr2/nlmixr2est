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
    f <- .nlmixr(one.cmt, theo_sd, "focei", foceiControl(maxOuterIterations = 0L))

    expect_null(f$foceiControl$etaMat)

    f2 <- .nlmixr(f, est = "focei", control = foceiControl(outerOpt = "bobyqa", maxOuterIterations = 0L))
    expect_true(inherits(f2$foceiControl$etaMat, "matrix"))

    f2 <- .nlmixr(f, "focei", foceiControl(outerOpt = "bobyqa", maxOuterIterations = 0L))
    expect_true(inherits(f2$foceiControl$etaMat, "matrix"))

    f3 <- .nlmixr(f, est = "focei", control = foceiControl(outerOpt = "bobyqa", etaMat = f, maxOuterIterations = 0L))
    expect_true(inherits(f3$foceiControl$etaMat, "matrix"))

    f4 <- .nlmixr(f, est = "focei", control = foceiControl(outerOpt = "bobyqa", etaMat = NA, maxOuterIterations = 0L))
    expect_true(is.na(f4$foceiControl$etaMat))

    f4 <- .nlmixr(f, "focei", foceiControl(outerOpt = "bobyqa", etaMat = NA, maxOuterIterations = 0L))

    expect_true(is.na(f4$foceiControl$etaMat))

    f4 <- .nlmixr(f, foceiControl(outerOpt = "bobyqa", etaMat = NA, maxOuterIterations = 0L))

    expect_true(is.na(f4$foceiControl$etaMat))
  })

  test_that("etaMat drops the mixture bookkeeping columns", {
    ## $eta carries ID, and for a mixture fit nmObjGet.ranef merges in a mixnum
    ## column.  Neither is an eta.  Letting either reach an etaMat trips
    ## foceiSetup_'s column check ("The etaMat must have the same number of ETAs
    ## (cols) as the model"), which is what broke every mixture refit -- $cov,
    ## addCwres, .setOfvFo and nlmixr2(fit, ...) all round-trip fit$etaMat.
    .eta <- data.frame(ID = 1:3, eta.ka = c(0.1, -0.2, 0.3), eta.cl = c(-0.1, 0.2, 0.0), mixnum = c(1L, 2L, 1L))
    expect_equal(colnames(.nmDropNonEtaCols(.eta)), c("eta.ka", "eta.cl"))
    ## the accessor itself, with and without the mixture column
    expect_equal(colnames(nmObjGet.etaMat(list(list(eta = .eta, iov = NULL)))), c("eta.ka", "eta.cl"))
    expect_equal(colnames(nmObjGet.etaMat(list(list(eta = .eta[, 1:3], iov = NULL)))), c("eta.ka", "eta.cl"))
    ## MIXEST is the same bookkeeping under its focei spelling
    .eta2 <- .eta
    names(.eta2)[4] <- "MIXEST"
    expect_equal(colnames(.nmDropNonEtaCols(.eta2)), c("eta.ka", "eta.cl"))
    ## an ordinary eta whose name merely contains "mix" is NOT dropped
    .eta3 <- data.frame(ID = 1:2, eta.mixup = c(0.1, 0.2))
    expect_equal(colnames(.nmDropNonEtaCols(.eta3)), "eta.mixup")
  })

  test_that("etaMat holds the occasion etas on the model's scale", {
    ## $iov reports the occasion etas times their standard deviation; the
    ## expanded model's occasion etas have unit variance
    .iovMod <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; iov.cl ~ 0.04 | occ })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl + iov.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd) })
    }
    .d <- nlmixr2data::theo_md
    .d$occ <- 1L + (.d$TIME >= 144)
    f <- .nlmixr(.iovMod, .d, "focei", foceiControl(print = 0L, maxOuterIterations = 0L, covMethod = ""))
    expect_equal(colnames(f$etaMat), c("eta.ka", "eta.cl", "eta.v", "rx.iov.cl.1", "rx.iov.cl.2"))
    ## held fixed, the fit's own etas reproduce its objective
    .ctl <- foceiControl(
      print = 0L,
      maxOuterIterations = 0L,
      maxInnerIterations = 0L,
      covMethod = "",
      etaMat = f$etaMat
    )
    expect_equal(.nlmixr(f$ui, .d, "focei", .ctl)$objf, f$objf, tolerance = 1e-8)
    ## imp replaces etaObf with a FOCEi recompute, and SAEM's two-level etas are
    ## natural-scale; $etaMat follows $eta and $iov for both
    for (f in list(
      .nlmixr(.iovMod, .d, "imp", impControl(print = 0L, nIter = 2L, covMethod = "")),
      .nlmixr(.iovMod, .d, "saem", saemControl(print = 0L, nBurn = 10L, nEm = 10L, covMethod = ""))
    )) {
      .sd <- sqrt(f$ui$omega$occ[1, 1])
      expect_equal(unname(f$etaMat[, 1:3]), unname(as.matrix(f$eta[, -1])))
      expect_equal(as.vector(t(f$etaMat[, 4:5])) * .sd, f$iov$occ$iov.cl)
    }
    ## a correlated block is expanded occasion by occasion, unscaled
    .corMod <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; iov.cl + iov.v ~ c(0.1, 0.03, 0.2) | occ })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl + iov.cl); v <- exp(tv + eta.v + iov.v)
        linCmt() ~ add(add.sd) })
    }
    f <- .nlmixr(.corMod, .d, "focei", foceiControl(print = 0L, maxOuterIterations = 0L, covMethod = ""))
    expect_equal(colnames(f$etaMat)[4:7], c("rx.iov.cl.1", "rx.iov.v.1", "rx.iov.cl.2", "rx.iov.v.2"))
    .ctl <- foceiControl(
      print = 0L,
      maxOuterIterations = 0L,
      maxInnerIterations = 0L,
      covMethod = "",
      etaMat = f$etaMat
    )
    expect_equal(.nlmixr(f$ui, .d, "focei", .ctl)$objf, f$objf, tolerance = 1e-8)
  })

  test_that(".foceiGradDirect() refits with an etaMat of etas only", {
    .mixMod <- function() {
      ini({ tka <- 0.45; tcl1 <- log(c(0, 2.7, 100)); tcl2 <- log(c(0, 0.1, 120)); tv <- 3.45; p1 <- 0.3
        add.sd <- 0.7; eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1 })
      model({ ka <- exp(tka + eta.ka); cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl)); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd) })
    }
    f <- .nlmixr(
      .mixMod,
      nlmixr2data::theo_sd,
      "focei",
      foceiControl(print = 0L, maxOuterIterations = 0L, covMethod = "")
    )
    expect_true("mixnum" %in% names(f$eta))
    .acc <- new.env(parent = emptyenv())
    local_mocked_bindings(nlmixr2 = function(object, data, est, control, ...) {
      .acc$etaMat <- control$etaMat
      NULL
    })
    .foceiGradDirect(f)
    expect_equal(colnames(.acc$etaMat), c("eta.ka", "eta.cl", "eta.v"))
  })
})
