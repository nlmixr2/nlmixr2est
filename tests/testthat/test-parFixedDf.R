nmTest({
  test_that("SD is consistent between parFixedDf and parFixed", {
    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- exp(0.45) # Log Ka
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
        ka <- tka + eta.ka
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    fit <- nlmixr(one.cmt, theo_sd, est="focei",
                  control=list(print=0))


    expect_equal(fit$parFixedDf[["BSV(CV% or SD)"]][1],
                 sqrt(fit$omega[1, 1]))

    # parFixed prints the SD to the foceiControl() default `sigdig`
    # significant figures (now 4); parFixedDf keeps full precision.
    expect_equal(as.numeric(fit$parFixed[["BSV(CV% or SD)"]][1]),
                 signif(sqrt(fit$omega[1, 1]), foceiControl()$sigdig))
  })
})
