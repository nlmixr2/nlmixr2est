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
  
  of1 <-f$objf
  etaMat1 <- as.matrix(f$eta[,-1])
  theta1 <- f$theta
  omega1 <- f$omega

  etaO1 <- f$etaObf
  
  f <- nlmixr(one.cmt, theo_sd, "foce")
  expect_true("CWRES" %in% names(f))

  of2 <- f$objf

  etaMat2 <- as.matrix(f$eta[,-1])
  theta2 <- f$theta
  omega2 <- f$omega

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
      # Use same form as stan math
      # https://github.com/stan-dev/math/blob/231cbd0be4280eb95cbb367068573830a7feae50/stan/math/prim/prob/normal_lpdf.hpp#L76-L80
      invSigma <- 1/add.sd
      yScaled <- (DV-cp)*invSigma
      yScaledSq <- yScaled*yScaled
      ll <- -0.5*yScaledSq - 0.5*log(2*pi) - log(add.sd)
      ll(err) ~ ll
    })
  }

  

  ## one.cmt.ll %>%
  ##   ini(theta1) %>%
  ##   ini(omega1) ->
  ##   one.cmt.ll

  ## foceiControl(etaMat=etaMat1, maxInnerIterations=0)
  
  f <- try(nlmixr(one.cmt.ll, theo_sd, "focei"))
  
  expect_false(inherits(f, "try-error"))
  expect_equal(-2*f$ll, f$IPRED)
  expect_false("CWRES" %in% names(f))

  expect_equal(f$objf, of1)

  f <- try(nlmixr(one.cmt.ll, theo_sd, "foce"))
  
  expect_false(inherits(f, "try-error"))
  expect_equal(-2*f$ll, f$IPRED)
  expect_false("CWRES" %in% names(f))

  expect_equal(f$objf, of2)
  
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
      ll(err) <- -log(add.sd*sqrt(2*pi))-0.5*((DV-cp)/add.sd)^2
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
      error ~ add(pdadd.err) + dnorm()
    })
  }

  f <- nlmixr(pk.turnover.emax3, nlmixr2data::warfarin, "focei")

  pk.turnover.emax3.ll <- function() {
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
      ll(pca) ~ -log(pdadd.err*sqrt(2*pi))-0.5*((DV-effect)/pdadd.err)^2
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

  pk.turnover.emax3.noeta <- pk.turnover.emax3 %>%
    model(cl <-exp(tcl)) %>%
    model(v <- exp(tv)) %>%
    model(emax <- expit(temax)) %>%
    model(kout <- exp(tkout)) %>%
    model(e0 <- exp(te0)) %>%
    ini(pkadd.err=fix(0.001))

  f <- nlmixr(pk.turnover.emax3.noeta, nlmixr2data::warfarin, "focei")


  f2 <- f %>% dplyr::filter(CMT == "pca")
  f1 <- f %>% dplyr::filter(CMT != "pca")


  for (i in c("IRES", "IWRES")) {
    expect_false(any(is.na(f1[[i]])))
    expect_true(all(is.na(f2[[i]])))
  }
  
})
