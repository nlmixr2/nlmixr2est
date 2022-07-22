test_that("test focei llik", {
  # dnorm() works

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
      linCmt() ~ add(add.sd) + dnorm()
    })
  }

  f <- nlmixr(one.cmt, theo_sd, "focei")
  expect_true("CWRES" %in% names(f))
  
  f <- nlmixr(one.cmt, theo_sd, "foce")
  expect_true("CWRES" %in% names(f))

  expect_error(nlmixr(one.cmt, theo_sd, "fo"))

  one.cmt.ll <- function() {
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
      cp <- linCmt()
      ll <- -log(add.sd*sqrt(2*pi))-0.5*((DV-cp)/add.sd)^2
      n2ll(err) ~ -2*ll
    })
  }

  f <- try(nlmixr(one.cmt.ll, theo_sd, "focei"))
  
  expect_false(inherits(f, "try-error"))
  expect_equal(-2*f$ll, f$IPRED)
  expect_false("CWRES" %in% names(f))

  f <- try(nlmixr(one.cmt.ll, theo_sd, "foce"))
  
  expect_false(inherits(f, "try-error"))
  expect_equal(-2*f$ll, f$IPRED)
  expect_false("CWRES" %in% names(f))
  
  # no etas test
  one.cmt.noeta <- function() {
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
      linCmt() ~ add(add.sd) + dnorm()
    })
  }

  f <- nlmixr(one.cmt.noeta, theo_sd, "focei")

  one.cmt.ll.noeta <- function() {
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
      cp <- linCmt()
      ll <- -log(add.sd*sqrt(2*pi))-0.5*((DV-cp)/add.sd)^2
      n2ll(err) ~ -2*ll
    })
  }

  f <- nlmixr(one.cmt.ll.noeta, theo_sd, "focei")


  pk.turnover.emax3 <- function() {
    ini({
      tktr <- log(1)
      tka <- log(1)
      tcl <- log(0.1)
      tv <- log(10)
      ##
      eta.cl ~ 2
      eta.v ~ 1
      prop.err <- 0.1
      pkadd.err <- 0.1
      ##
      temax <- logit(0.8)
      tec50 <- log(0.5)
      tkout <- log(0.05)
      te0 <- log(100)
      ##
      eta.emax ~ .5
      eta.kout ~ .5
      eta.e0 ~ .5
      ##
      pdadd.err <- 10
    })
    model({
      ktr <- exp(tktr)
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      emax = expit(temax+eta.emax)
      ec50 =  exp(tec50)
      kout = exp(tkout + eta.kout)
      e0 = exp(te0 + eta.e0)
      ##
      DCP = center/v
      PD=1-emax*DCP/(ec50+DCP)
      ##
      effect(0) = e0
      kin = e0*kout
      ##
      d/dt(depot) = -ktr * depot
      d/dt(gut) =  ktr * depot -ka * gut
      d/dt(center) =  ka * gut - cl / v * center
      d/dt(effect) = kin*PD -kout*effect
      ##
      cp = center / v
      cp ~ prop(prop.err) + add(pkadd.err)
      lle <- -2*(-log(pdadd.err*sqrt(2*pi))-0.5*((DV-effect)/pdadd.err)^2)
      n2ll(pca) ~ lle
    })
  }
    
  f <- nlmixr(pk.turnover.emax3, nlmixr2data::warfarin, "focei")

  f1 <- f %>% dplyr::filter(CMT != "pca")
  f2 <- f %>% dplyr::filter(CMT == "pca")

  for (i in c("RES", "WRES", "IRES", "IWRES", "RES", "WRES",
              "IRES", "IWRES", "CPRED", "CRES", "CWRES")) {
    expect_false(any(is.na(f1[[i]])))
    expect_true(all(is.na(f2[[i]])))
  }

  f <- nlmixr(pk.turnover.emax3, nlmixr2data::warfarin, "foce")

  f2 <- f %>% dplyr::filter(CMT == "pca")
  f1 <- f %>% dplyr::filter(CMT != "pca")


  for (i in c("RES", "WRES", "IRES", "IWRES", "RES", "WRES",
              "IRES", "IWRES", "CPRED", "CRES", "CWRES")) {
    expect_false(any(is.na(f1[[i]])))
    expect_true(all(is.na(f2[[i]])))
  }
  

  

  
})
