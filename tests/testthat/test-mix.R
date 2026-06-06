nmTest({

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
