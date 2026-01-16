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

      f <- suppressMessages(suppressWarnings(nlmixr(One.SD.ODE, dat, "posthoc")))

      expect_s3_class(f, "nlmixr2.posthoc")

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
    fit <- .nlmixr(one.compartment, theo_sd, est = "focei", control=list(print=0))
    expect_s3_class(fit, "nlmixr2FitCore")
  })
})
