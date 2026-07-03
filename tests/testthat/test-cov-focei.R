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

  test_that("covMethod r/s/r,s SEs match NONMEM on theo_sd (issues #666, #694)", {
    # Validated against NONMEM 7 (FOCEI, ADVAN2 TRANS2, additive error) on the
    # identical theo_sd fit: $COV MATRIX=R ("r" = R^-1), MATRIX=S ("s" = S^-1),
    # and the default sandwich R^-1 S R^-1 ("r,s").  Estimates and OFV land on the
    # same optimum; the SEs below are the NONMEM reference on nlmixr2's native
    # scale (theta log-scale; Omega on the variance scale; add.sd as an SD).
    #
    # This guards two things: (#666) the covariance SCALING -- a spurious 2x or
    # sqrt(2) factor would miss NONMEM by ~40%; and (#694) the residual (sigma)
    # and Omega SEs that the cov now reports.  Note the "s" (cross-product / OPG)
    # theta SEs are LARGE (e.g. tcl ~0.43) -- this is the genuine inflation of the
    # OPG estimator once variance components enter the score matrix, and NONMEM's
    # MATRIX=S reproduces it to <1%.  So "s" must match the (inflated) reference,
    # not be "small".  "r" and the sandwich are the robust methods.
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
    .nm <- list(
      "r"   = c(tka = 0.1886, tcl = 0.0836, tv = 0.0466, add.sd = 0.0504,
                om.tka = 0.1881, om.tcl = 0.0346, om.tv = 0.0111),
      "s"   = c(tka = 0.2962, tcl = 0.4270, tv = 0.2233, add.sd = 0.0404,
                om.tka = 0.2732, om.tcl = 0.0985, om.tv = 0.0298),
      "r,s" = c(tka = 0.1908, tcl = 0.0739, tv = 0.0428, add.sd = 0.0998,
                om.tka = 0.2258, om.tcl = 0.0299, om.tv = 0.0084))
    for (.m in names(.nm)) {
      # covFull=TRUE so $cov spans theta + sigma + Omega (the reference SEs below
      # include the residual and Omega blocks); the default covFull=FALSE keeps the
      # original theta-only $cov.
      .f <- .nlmixr(one.cmt, theo_sd, "focei",
                    foceiControl(covMethod = .m, sigdig = 4, print = 0, covFull = TRUE))
      .se <- sqrt(pmax(diag(.f$cov), 0))
      .ref <- .nm[[.m]]
      expect_true(all(names(.ref) %in% names(.se)),
                  label = paste0("covMethod=", .m, ": cov has theta + sigma + Omega names"))
      .rel <- abs(.se[names(.ref)] / .ref - 1)
      expect_true(all(.rel < 0.15),
                  label = paste0("covMethod=", .m, " SEs vs NONMEM (max rel err ",
                                 round(max(.rel), 3), ")"))
    }
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
