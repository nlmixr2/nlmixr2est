nmTest({
  test_that("input ui doesn't change", {

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

    ui <- .nlmixr(one.cmt)

    est0 <- ui$iniDf$est

    fit <- .nlmixr(ui, theo_sd, est="focei", control=list(print=0))

    expect_equal(est0, ui$iniDf$est)

    # Also check that there are no sensitives from the linCmt()
    expect_false(any(names(fit) == "rx__sens_central_BY_p1"))
  })
})
