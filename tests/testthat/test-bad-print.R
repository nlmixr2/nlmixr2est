nmTest({
  test_that("odd message does not work", {

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
    run <- function(){
      call <- '
f <- suppressMessages(nlmixr2(one.cmt, theo_sd, "focei",
  control=foceiControl(print=0, outerOpt="lbfgsb3c"),
    table = list(cwres = TRUE, npde = TRUE)))
'
      eval(parse(text=call))
      return(f)
    }

    tf <- tempfile()

    fit <- run()
    withr::with_output_sink(tf,
                            expect_error(print(fit), NA))
    unlink(tf)


  })
})
