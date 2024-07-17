nmTest({
  test_that("funny one-compartment model", {

    one.cmt <- function() {
      ini({
        tka <- 0.45 ; label("Ka")
        tcl <- log(c(0, 2.7, 100)) ; label("Log Cl")
        tv <- 3.45; label("log V")
        cl.wt <- 0
        v.wt <- 0
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + WT * cl.wt)
        v <- exp(tv + eta.v) + WT ^ 2 * v.wt
        linCmt() ~ add(add.sd)
      })
    }

    skip_if_not(rxode2::.linCmtSensB())


    expect_error(nlmixr(one.cmt, nlmixr2data::theo_sd, "saem"), NA)

    expect_error(nlmixr(one.cmt, nlmixr2data::theo_sd, "focei"), NA)

  })
})
