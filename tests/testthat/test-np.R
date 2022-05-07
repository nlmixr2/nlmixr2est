nmTest({

  test_that("Test naive pooled", {

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }

    expect_error(suppressMessages(suppressWarnings(nlmixr(one.cmt, nlmixr2data::theo_sd, "foce"))), NA)
    expect_error(suppressMessages(suppressWarnings(nlmixr(one.cmt, nlmixr2data::theo_sd, "focei"))), NA)

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        eta.ka ~ 0
        eta.cl ~ 0
        eta.v ~ 0
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    expect_error(nlmixr(one.cmt, nlmixr2data::theo_sd, "foce"), NA)

    expect_error(nlmixr(one.cmt, nlmixr2data::theo_sd, "focei"), NA)
    expect_error(nlmixr(one.cmt, nlmixr2data::theo_sd, "saem"))

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        eta.ka ~ 0
        eta.cl ~ 0.1
        eta.v ~ 0
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    expect_error(nlmixr(one.cmt, nlmixr2data::theo_sd, "foce"), NA)
    expect_error(nlmixr(one.cmt, nlmixr2data::theo_sd, "focei"), NA)
    expect_error(nlmixr(one.cmt, nlmixr2data::theo_sd, "saem"), NA)

  })

})
