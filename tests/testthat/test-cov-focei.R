nmTest({

  test_that("cov-focei", {

    dat <-
      warfarin |>
      dplyr::filter(dvid == "cp")

    #### doesn't work with FOCEI ### doesn't work with SAEM with CRWES=True but does iwth CWRES=FALSE
    One.SD.ODE <- function() {
      ini({
        # Where initial conditions/variables are specified
        lcl <- log(0.135) # log Cl (L/h)
        lv <- log(8) # log V (L)
        lmtt <- log(1.1) # log MTx

        prop.err <- 0.15 # proportional error (SD/mean)
        add.err <- 0.6 # additive error (mg/L)
        eta.cl + eta.v + eta.mtt ~ c(
          0.1,
          0.001, 0.1,
          0.001, 0.001, 0.1
        )
      })
      model({
        # Where the model is specified
        cl <- exp(lcl + eta.cl)
        v <- exp(lv + eta.v)
        mtt <- exp(lmtt + eta.mtt)
        ktr <- 6 / mtt

        ## ODE example
        d / dt(depot) <- -ktr * depot
        d / dt(central) <- ktr * trans5 - (cl / v) * central
        d / dt(trans1) <- ktr * depot - ktr * trans1
        d / dt(trans2) <- ktr * trans1 - ktr * trans2
        d / dt(trans3) <- ktr * trans2 - ktr * trans3
        d / dt(trans4) <- ktr * trans3 - ktr * trans4
        d / dt(trans5) <- ktr * trans4 - ktr * trans5

        ## Concentration is calculated
        cp <- central / v
        ## And is assumed to follow proportional and additive error
        cp ~ prop(prop.err) + add(add.err)
      })
    }

    f <- .nlmixr(One.SD.ODE, dat, "posthoc")

    expect_s3_class(f, "nlmixr2.posthoc")
  })

  test_that("covMethod r and s SEs are not inflated vs the sandwich (issue #666)", {
    one.cmt <- function() {
      ini({
        tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    fit_r  <- .nlmixr(one.cmt, theo_sd, "focei", foceiControl(covMethod = "r",   print = 0))
    fit_s  <- .nlmixr(one.cmt, theo_sd, "focei", foceiControl(covMethod = "s",   print = 0))
    fit_rs <- .nlmixr(one.cmt, theo_sd, "focei", foceiControl(covMethod = "r,s", print = 0))

    p <- c("tka", "tcl", "tv")
    se_r  <- fit_r$parFixedDf[p,  "SE"]
    se_s  <- fit_s$parFixedDf[p,  "SE"]
    se_rs <- fit_rs$parFixedDf[p, "SE"]

    # Before the #666 fix: SE_r was ~sqrt(2)*SE_rs (covR was 2*Rinv) and SE_s was ~2x that again
    # (covS was 4*Sinv).  covMethod="r" (observed information) and the "r,s" sandwich obey the
    # information equality, so r/rs stays ~1 and a return of the sqrt(2) scaling would push it >1.41.
    expect_true(all(se_r  / se_rs < 1.3), label = "covMethod='r' SE not inflated vs sandwich")
    # covMethod="s" is the OPG (score cross-product) estimator.  On this 12-subject dataset it
    # legitimately disagrees with the Hessian in its off-diagonals (S corr(tcl,tv) ~0.85 vs R's
    # ~0.12), so s/rs runs ~2 for cl/v -- finite-sample OPG behaviour that collapses toward 1 on
    # larger / well-specified data, NOT a scaling error.  This leg guards against the old constant
    # factor (covS = 4*Sinv), which would DOUBLE these ratios to ~4.5.
    expect_true(all(se_s  / se_rs < 3), label = "covMethod='s' SE not inflated by the old 2x constant factor")
  })

  test_that("covariance with many omegas fixed will not crash focei", {
    one.compartment <- function() {
      ini({
        tka <- fix(0.45) # Log Ka
        tcl <- 1 # Log Cl
        tv <- fix(3.45)    # Log V
        eta.ka ~ fix(0.6)
        eta.cl ~ fix(0.3)
        eta.v ~ fix(0.1)
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.sd)
      })
    }
    fit <- .nlmixr(one.compartment, theo_sd, est = "focei",
                   control=foceiControl(print=0, maxOuterIterations=0L))
    expect_s3_class(fit, "nlmixr2FitCore")
  })
})
