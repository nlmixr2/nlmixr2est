nmTest({
  test_that("vpcSim retains column information", {

    PKPDdata <- nlmixr2data::warfarin
    PKPDdata$dvid <- ifelse(PKPDdata$dvid == "cp", "central", "effect")
    PKPDdata <- PKPDdata[, names(PKPDdata) != "DVID"]

    KA1tr1_PDimmemax1 <- function() {
      ini({
        ## PK
        tktr <- log(1)   # log ktr (/h)
        tcl  <- log(0.1) # log CL (L/h)
        tv   <- log(8)   # log Vc (L)

        eta.ktr ~ 1
        eta.cl ~ 0.1
        eta.v ~ 0.1
        eps.pkprop <- 0.1
        eps.pkadd <- 0.4

        ## PD
        tc50  <- log(1)    #log ec50 (mg/L)
        te0   <- log(100)  #log e0

        eta.c50  ~ .5
        eta.e0 ~ .1

        eps.pdadd <- 100

      })
      model({
        ktr <- exp(tktr + eta.ktr)
        cl  <- exp(tcl + eta.cl)
        v   <- exp(tv + eta.v)

        c50  = exp(tc50 + eta.c50)
        e0   = exp(te0 + eta.e0)

        cp           =  central/v
        d/dt(depot)  = -ktr * depot
        d/dt(central)=  ktr * trans - cl * cp
        d/dt(trans)  =  ktr * depot - ktr * trans
        effect       =  e0 * (1 - cp/(c50 + cp))

        cp ~ prop(eps.pkprop) + add(eps.pkadd) | central
        effect ~ add(eps.pdadd) | effect
      })
    }

    fitKA1tr1_PDimmemax1_F <-
      nlmixr(KA1tr1_PDimmemax1,
             PKPDdata,
             est = "focei",
             foceiControl(print = 5))

    f <- vpcSim(fitKA1tr1_PDimmemax1_F)

    expect_true(inherits(f$CMT, "factor"))

    tmp <- vpcSimExpand(fitKA1tr1_PDimmemax1_F, f, "dvid")

    expect_true(inherits(tmp$dvid, "factor"))

    tmp2 <- tmp

    tmp2$CMT <- paste(tmp2$CMT)
    tmp2$dvid <- paste(tmp2$dvid)

    tmp2 <- vpcNameDataCmts(fitKA1tr1_PDimmemax1_F, tmp2)

    expect_equal(tmp$CMT, tmp2$CMT)

    expect_equal(tmp$dvid, tmp2$dvid)

  })

  test_that("vpcSim works with etas that are set to zero (#341)", {
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- log(c(0, 2.7, 100))
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
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

    skip_if_not(rxode2parse::.linCmtSens())

    suppressMessages(
      fit <- nlmixr(one.cmt, theo_sd, est="focei", control = foceiControl(print = 0, eval.max = 1))
    )
    expect_s3_class(vpcSim(fit, pred=TRUE), "data.frame")
  })
})
