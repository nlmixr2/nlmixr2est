nmTest({

  test_that("vpcSim works for IOV models (simInfo omega includes per-occasion etas)", {
    one.cmt.iov <- function() {
      ini({
        tka <- 0.45
        tcl <- log(c(0, 2.7, 100))
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        iov.cl ~ 0.1 | occ
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + iov.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    theo_iov <- nlmixr2data::theo_md
    theo_iov$occ <- 1L
    theo_iov$occ[theo_iov$TIME >= 144] <- 2L

    fit <- .nlmixr(one.cmt.iov, theo_iov, est = "focei",
                   control = foceiControlFast)

    # simInfo omega must include per-occasion IOV etas
    .si <- fit$simInfo
    .omegaNames <- dimnames(.si$omega)[[1]]
    expect_true(any(grepl("^rx\\.iov\\.cl\\.", .omegaNames)),
                info = "simInfo omega should contain per-occasion IOV etas (rx.iov.cl.*)")

    # Two occasions → two per-occasion etas
    .iovCols <- grep("^rx\\.iov\\.cl\\.", .omegaNames, value = TRUE)
    expect_equal(length(.iovCols), 2L,
                 info = "Two occasions should produce two per-occasion IOV etas in omega")

    # params must include the IOV theta (iov.cl)
    expect_true("iov.cl" %in% names(.si$params),
                info = "simInfo params should include the IOV theta parameter")

    # vpcSim should run without error and return a data frame
    vpc <- vpcSim(fit, n = 50, seed = 42)
    expect_s3_class(vpc, "data.frame")
    expect_true(nrow(vpc) > 0)

    # pred=TRUE also works
    vpc2 <- vpcSim(fit, n = 50, seed = 42, pred = TRUE)
    expect_s3_class(vpc2, "data.frame")
    expect_true("pred" %in% names(vpc2))
  })

  test_that("vpcSim IOV: iovXform is persisted on fit and used in simulation", {
    one.cmt.iov <- function() {
      ini({
        tka <- 0.45
        tcl <- log(c(0, 2.7, 100))
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        iov.cl ~ 0.1 | occ
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + iov.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    theo_iov <- nlmixr2data::theo_md
    theo_iov$occ <- 1L
    theo_iov$occ[theo_iov$TIME >= 144] <- 2L

    .ctl <- foceiControlFast
    .ctl$iovXform <- "logvar"
    fit <- .nlmixr(one.cmt.iov, theo_iov, est = "focei", control = .ctl)

    # iovXform should be stored on the fit environment
    expect_true(exists("iovXform", fit$env),
                info = "iovXform should be stored in fit$env after estimation")
    expect_equal(fit$env$iovXform, "logvar")

    # Simulation should still work with non-default xform
    vpc <- vpcSim(fit, n = 50, seed = 42)
    expect_s3_class(vpc, "data.frame")
    expect_true(nrow(vpc) > 0)
  })

  test_that("vpcSim non-IOV models unaffected by IOV changes", {
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- log(c(0, 2.7, 100))
        tv <- 3.45
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

    fit <- .nlmixr(one.cmt, theo_sd, est = "focei",
                   control = foceiControl(print = 0, eval.max = 1))

    # No iovXform on non-IOV fit
    expect_false(exists("iovXform", fit$env),
                 info = "Non-IOV fit should not have iovXform in env")

    # omega should only have BSV etas (3)
    .si <- fit$simInfo
    expect_equal(nrow(.si$omega), 3L,
                 info = "Non-IOV model should have exactly 3 etas in omega")

    vpc <- vpcSim(fit, n = 50, seed = 42)
    expect_s3_class(vpc, "data.frame")
    expect_true(nrow(vpc) > 0)
  })

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
      .nlmixr(
        KA1tr1_PDimmemax1,
        PKPDdata,
        est = "focei",
        foceiControl(print = 5)
      )

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

    fit <- .nlmixr(one.cmt, theo_sd, est="focei", control = foceiControl(print = 0, eval.max = 1))
    expect_s3_class(vpcSim(fit, pred=TRUE), "data.frame")
  })
})
