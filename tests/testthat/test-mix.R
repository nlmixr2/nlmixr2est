nmTest({

  test_that("test mixture models -- focei fit", {
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl1 <- log(c(0, 2.7, 100))
        tcl2 <- log(c(0, 0.1, 120))
        tv <- 3.45
        p1 <- 0.3
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    withr::with_seed(42, {
      fit <- nlmixr2(one.cmt, nlmixr2data::theo_sd, "focei",
                     control = foceiControl(print = 0, maxOuterIterations = 5))
    })

    # ranef should have ID + ETA columns only (no MIXEST)
    expect_true(all(c("ID", "eta.ka", "eta.cl", "eta.v") %in% names(fit$ranef)))
    expect_false("MIXEST" %in% names(fit$ranef))
    expect_equal(nrow(fit$ranef), length(unique(nlmixr2data::theo_sd$ID)))

    # mixNum: one row per subject
    mn <- fit$mixNum
    expect_true(is.data.frame(mn))
    expect_true(all(c("ID", "MIXNUM") %in% names(mn)))
    expect_equal(nrow(mn), length(unique(nlmixr2data::theo_sd$ID)))
    expect_true(all(mn$MIXNUM %in% c(1L, 2L)))

    # mixList: one component per mixture
    ml <- fit$mixList
    expect_true(is.list(ml))
    expect_equal(length(ml), 2L)
    expect_true(all(c("mix1", "mix2") %in% names(ml)))
    for (k in seq_along(ml)) {
      expect_true(all(c("ID", "eta.ka", "eta.cl", "eta.v", "prob") %in% names(ml[[k]])))
      expect_equal(nrow(ml[[k]]), length(unique(nlmixr2data::theo_sd$ID)))
      expect_true(all(ml[[k]]$prob >= 0 & ml[[k]]$prob <= 1))
    }

    # Residuals/table step completed (CWRES exists)
    expect_true("CWRES" %in% names(fit))
  })

  test_that("test mixture models -- ui components", {

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl1 <- log(c(0, 2.7, 100)) # Log Cl
        tcl2 <- log(c(0, 0.1, 120)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        p1 <- 0.3
        ## the label("Label name") works with all models
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v <- exp(tv + eta.v)
        me <- mixest
        mn <- mixnum
        mu <- mixunif
        linCmt() ~ add(add.sd)
      })
    }

    ui <- rxode2::rxode2(one.cmt())

    expect_equal(ui$thetaIniMix,
                 c(tka = 0.45, tcl1 = 0.993251773010283, tcl2 = -2.30258509299405,
                   tv = 3.45, p1 = -0.847297860387204, add.sd = 0.7))

    expect_equal(ui$thetaMixIndex, 5L)


    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl1 <- log(c(0, 2.7, 100)) # Log Cl
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
        cl <- exp(tcl1 + eta.cl)
        v <- exp(tv + eta.v)
        me <- mixest
        mn <- mixnum
        mu <- mixunif
        linCmt() ~ add(add.sd)
      })
    }

    ui <- rxode2::rxode2(one.cmt())

    expect_equal(ui$thetaIniMix,
                 c(tka = 0.45, tcl1 = log(2.7), tv = 3.45, add.sd = 0.7))

    expect_equal(ui$thetaMixIndex, integer(0))

  })

})
