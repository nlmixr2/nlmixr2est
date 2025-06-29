nmTest({
  test_that("ipred=1 will not occur", {

    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        Nominal = Nominal
        linCmt() ~ add(add.sd)
      })
    }
    f <- nlmixr(one.cmt)

    theo_sd2 <- nlmixr2data::theo_sd
    theo_sd2$Nominal <- 1
    fit <- nlmixr(one.cmt, theo_sd2, est="nlm", control=list(print=0))

    expect_true(all(fit$IPRED != 1))

  })
})
