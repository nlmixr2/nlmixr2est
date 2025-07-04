nmTest({
  test_that("rxode2 model properties work with fits", {

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Ka
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

    f <- .nlmixr(one.cmt)

    fit <- .nlmixr(one.cmt, theo_sd, est="saem", control = saemControlFast)

    expect_equal(rxModelVars(one.cmt), rxModelVars(f))

    expect_equal(rxModelVars(f), rxModelVars(fit))
  })
})
