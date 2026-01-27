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

    expect_error(.nlmixr(one.cmt, nlmixr2data::theo_sd, "foce", list(print=0)), NA)
    expect_error(.nlmixr(one.cmt, nlmixr2data::theo_sd, "focei", list(print=0)), NA)

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

    expect_error(.nlmixr(one.cmt, nlmixr2data::theo_sd, "foce", list(print=0)), NA)

    expect_error(.nlmixr(one.cmt, nlmixr2data::theo_sd, "focei", list(print=0)), NA)
    expect_error(
      .nlmixr(one.cmt, nlmixr2data::theo_sd, "saem", list(print=0)),
      regexp = "needs to be a mixed effect model for the estimation routine 'saem'",
      fixed = TRUE
    )

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

    expect_error(.nlmixr(one.cmt, nlmixr2data::theo_sd, "foce", list(print=0)), NA)
    expect_error(.nlmixr(one.cmt, nlmixr2data::theo_sd, "focei", list(print=0)), NA)
    expect_error(.nlmixr(one.cmt, nlmixr2data::theo_sd, "saem", control = saemControlFast), NA)
  })
})

test_that("parameters are updated in fit object", {

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

  fit <- suppressMessages(suppressWarnings(.nlmixr(one.cmt, nlmixr2data::theo_sd, "focei",
                                                  list(print=0))))

  f1 <- suppressMessages(one.cmt())

  expect_false(isTRUE(all.equal(f1$iniDf$est, fit$iniDf$est)))

})
