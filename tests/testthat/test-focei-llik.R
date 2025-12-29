if (rxode2hasLlik()) {
  nmTest({
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

      f <- .nlmixr(one.cmt, theo_sd, "focei")
      expect_true("CWRES" %in% names(f))

      of1     <- f$objf
      etaMat1 <- as.matrix(f$eta[,-1])
      theta1  <- f$theta
      omega1  <- f$omega
      etaO1   <- f$etaObf

      f <- .nlmixr(one.cmt, theo_sd, "foce")
      expect_true("CWRES" %in% names(f))

      of2 <- f$objf
      etaMat2 <- as.matrix(f$eta[,-1])
      theta2 <- f$theta
      omega2 <- f$omega

      expect_error(.nlmixr(one.cmt, theo_sd, "fo"))

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

      one.cmt.ll |>
        ini(theta1) |>
        ini(omega1) ->
        one.cmt.ll

      f <- try(.nlmixr(one.cmt.ll, theo_sd, "focei",
                      control=foceiControl(etaMat=etaMat1, maxInnerIterations=0,
                                           maxOuterIterations=0)))

      expect_false(inherits(f, "try-error"))
      expect_equal(f$ll, f$IPRED)
      expect_false("CWRES" %in% names(f))

      expect_equal(f$objf, of1)

      one.cmt.ll |>
        ini(theta2) |>
        ini(omega2) ->
        one.cmt.ll

      f <- try(.nlmixr(one.cmt.ll, theo_sd, "foce",
                      control=foceiControl(etaMat=etaMat2, maxInnerIterations=0,
                                           maxOuterIterations=0)))

      expect_false(inherits(f, "try-error"))
      expect_equal(f$ll, f$IPRED)
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

      f <- .nlmixr(one.cmt.noeta, theo_sd, "focei")

      of1 <-f$objf
      theta1 <- f$theta

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
          ll(err) ~ -log(add.sd) - 0.5*log(2*pi) - 0.5*((DV-cp)/add.sd)^2
        })
      }

      one.cmt.ll.noeta |>
        ini(theta1) ->
        one.cmt.ll.noeta

      f <- .nlmixr(one.cmt.ll.noeta, theo_sd, "focei",
                  control=foceiControl(maxOuterIterations=0))

      expect_equal(of1, f$objf)
    })


    pk.turnover.emax3.n1 <- function() {
      ini({
        tktr <- log(1)
        tka <- log(1)
        tcl <- log(0.1)
        tv <- log(10)
        ##
        eta.ktr ~ 1
        eta.ka ~ 1
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
        eta.ec50  ~ .5
        eta.kout ~ .5
        eta.e0 ~ .5
        ##
        pdadd.err <- 1
      })
      model({
        ktr <- exp(tktr + eta.ktr)
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        emax = expit(temax+eta.emax)
        ec50 =  exp(tec50 + eta.ec50)
        kout = exp(tkout + eta.kout)
        e0 = exp(te0 + eta.e0)
        ##
        DCP = center/v
        rPD=1-emax*DCP/(ec50+DCP)
        ##
        effect(0) = e0
        kin = e0*kout
        ##
        d/dt(depot) = -ktr * depot
        d/dt(gut) =  ktr * depot -ka * gut
        d/dt(center) =  ka * gut - cl / v * center
        d/dt(effect) = kin*rPD -kout*effect
        ##
        cp = center / v
        cp ~ prop(prop.err) + add(pkadd.err)
        effect ~ add(pdadd.err) + dnorm() | pca
      })
    }

    pk.turnover.emax3.n2 <- function() {
      ini({
        tktr <- log(1)
        tka <- log(1)
        tcl <- log(0.1)
        tv <- log(10)
        ##
        eta.ktr ~ 1
        eta.ka ~ 1
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
        eta.ec50  ~ .5
        eta.kout ~ .5
        eta.e0 ~ .5
        ##
        pdadd.err <- 1
      })
      model({
        ktr <- exp(tktr + eta.ktr)
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        emax = expit(temax+eta.emax)
        ec50 =  exp(tec50 + eta.ec50)
        kout = exp(tkout + eta.kout)
        e0 = exp(te0 + eta.e0)
        ##
        DCP = center/v
        rPD=1-emax*DCP/(ec50+DCP)
        ##
        effect(0) = e0
        kin = e0*kout
        ##
        d/dt(depot) = -ktr * depot
        d/dt(gut) =  ktr * depot -ka * gut
        d/dt(center) =  ka * gut - cl / v * center
        d/dt(effect) = kin*rPD -kout*effect
        ##
        cp = center / v
        cp ~ prop(prop.err) + add(pkadd.err) + dnorm()
        effect ~ add(pdadd.err) + dnorm() | pca
      })
    }

    f <- .nlmixr(pk.turnover.emax3.n1)
    f2 <- .nlmixr(pk.turnover.emax3.n2)

    expect_equal(f$foceiModel0, f2$foceiModel0)

    expect_equal(f$foceModel0, f2$foceModel0)

    f <- .nlmixr(pk.turnover.emax3.n1, nlmixr2data::warfarin, "focei",
                control=foceiControl(covMethod = "",
                                     maxOuterIterations=0))

    of1     <- f$objf
    etaMat1 <- as.matrix(f$eta[,-1])
    theta1  <- f$theta
    omega1  <- f$omega
    etaO1   <- f$etaObf


    pk.turnover.emax3.ll <- function() {
      ini({
        tktr <- log(1)
        tka <- log(1)
        tcl <- log(0.1)
        tv <- log(10)
        ##
        eta.ktr ~ 1
        eta.ka ~ 1
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
        eta.ec50  ~ .5
        eta.kout ~ .5
        eta.e0 ~ .5
        ##
        pdadd.err <- c(0, 10)
      })
      model({
        ktr <- exp(tktr + eta.ktr)
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        emax = expit(temax+eta.emax)
        ec50 =  exp(tec50 + eta.ec50)
        kout = exp(tkout + eta.kout)
        e0 = exp(te0 + eta.e0)
        ##
        DCP = center/v
        rPD=1-emax*DCP/(ec50+DCP)
        ##
        effect(0) = e0
        kin = e0*kout
        ##
        d/dt(depot) = -ktr * depot
        d/dt(gut) =  ktr * depot -ka * gut
        d/dt(center) =  ka * gut - cl / v * center
        d/dt(effect) = kin*rPD -kout*effect
        ##
        cp = center / v
        cp ~ prop(prop.err) + add(pkadd.err) + dnorm()
        ll(pca) ~ -log(pdadd.err) - 0.5*log(2*pi) - 0.5 * ((DV - effect) / pdadd.err)^2
      })
    }

    pk.turnover.emax3.ll |>
      ini(theta1) |>
      ini(omega1) ->
      pk.turnover.emax3.ll

    f2 <- .nlmixr(pk.turnover.emax3.ll, nlmixr2data::warfarin, "focei",
                 control=foceiControl(etaMat=etaMat1, maxInnerIterations=0,
                                      maxOuterIterations=0,
                                      optimHessType="forward"))

    test_that("same values for omega, theta and eta, forward", {
      expect_equal(f$omega, f2$omega)
      expect_equal(f$theta, f2$theta)
      expect_equal(f$eta, f2$eta)
    })

    f2 <- .nlmixr(pk.turnover.emax3.ll, nlmixr2data::warfarin, "focei",
                 control=foceiControl(etaMat=etaMat1, maxInnerIterations=0,
                                      maxOuterIterations=0,
                                      optimHessType="central"))

    test_that("same values for omega, theta and eta, central", {
      expect_equal(f$omega, f2$omega)
      expect_equal(f$theta, f2$theta)
      expect_equal(f$eta, f2$eta)
    })


    test_that("objective values are equal for mixed ll", {
      expect_equal(f2$objf, of1)
    })

    fll <- addNpde(f2)
    fnorm <- addNpde(f)

    f1 <- fll |> dplyr::filter(CMT != "pca")
    f1norm <- fnorm |> dplyr::filter(CMT != "pca")
    f2 <- fll |> dplyr::filter(CMT == "pca")

    for (i in c("RES", "WRES", "IRES", "IWRES", "WRES",
                "IWRES", "CPRED", "CRES", "CWRES", "PRED", "IPRED",
                "EPRED", "ERES", "NPDE", "NPD",
                "PDE", "PD")) {
      test_that(paste0("res: ", i), {
        expect_false(any(is.na(f1[[i]])))
        if (i %in% c("PRED", "IPRED")) {
          expect_false(any(is.na(f2[[i]])))
        } else {
          expect_true(all(is.na(f2[[i]])))
        }
        if (!(i %in% c("EPRED", "ERES", "NPDE", "NPD", "PDE", "PD"))) {
          expect_equal(f1[[i]], f1norm[[i]])
        }
      })
    }
  })
}
