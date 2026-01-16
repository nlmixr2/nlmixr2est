nmTest({
  mod <- function() {
    ini({
      tka <- 0.45 ; label("Log Ka")
      tcl <- 1 ; label("Log Cl")
      tv <- 3.45 ; label("Log V")
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
      wt.v2 <- -3
      sexf.cl <- 1.5
      wt.cl <- 3
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + wt.cl*logWt70 + sexf.cl*sexf)
      v <- exp(tv + eta.v + wt.v2*logWt70)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      ipre <- center / v
      ipre ~ add(add.sd)
    })
  }

  nid <- 12

  dat <- rxode2::rxWithSeed(102478, {
    cov <- data.frame(id=seq(nid), logWt70=log(rnorm(nid, 70, 10) / 70), sexf=round(runif(nid)))
    ev <- et() |>
      et(amt=320) |>
      et(c(0.25, 0.5, 1, 2, 4, 6, 8, 10, 12, 24)) |>
      et(id=1:nid) |>
      merge(cov)
    rxSolve(mod, ev, addDosing=TRUE, returnType="data.frame") |>
      dplyr::select(id, evid, cmt, amt, time, sim, logWt70, sexf) |>
      dplyr::rename(dv=sim)
  })

  test_that("test mu reference covariate in saem", {
    mod2 <- function() {
      ini({
        wt.cl <- 3
        sexf.cl <- 1.5
        wt.v2 <- -3
        tka <- 0.45 ; label("Log Ka")
        tcl <- 1 ; label("Log Cl")
        tv <- 3.45 ; label("Log V")
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + wt.cl*logWt70 + sexf.cl*sexf)
        v <- exp(tv + eta.v + wt.v2*logWt70)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        ipre <- center / v
        ipre ~ add(add.sd)
      })
    }

    f1 <- .nlmixr(mod)

    f2 <- .nlmixr(mod2)
    expect_equal(f1$saemParamsToEstimate, f2$saemParamsToEstimate)
    expect_equal(f1$saemInitTheta, f2$saemInitTheta)

    fit1 <- .nlmixr(mod, dat, "saem", control = saemControlFast)
    fit2 <- .nlmixr(mod2, dat, "saem", control = saemControlFast)

    for (n in names(fit1$theta)) {
      expect_equal(fit1$theta[n], fit2$theta[n])
    }

    expect_false(all(fit1$theta == fit2$theta))
  })

  test_that("test fixed parameters become non mu referenced covariates", {
    mod <- function() {
      ini({
        tka <- 0.45 ; label("Log Ka")
        tcl <- 1 ; label("Log Cl")
        tv <- 3.45 ; label("Log V")
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
        wt.v2 <- -3
        sexf.cl <- fix(1.5)
        wt.cl <- 3
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + wt.cl*logWt70 + sexf.cl*sexf)
        v <- exp(tv + eta.v + wt.v2*logWt70)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        ipre <- center / v
        ipre ~ add(add.sd)
      })
    }

    f <- .nlmixr(mod)

    expect_equal(f$saemParamsToEstimateCov, c("tka", "tcl", "tv", "sexf.cl"))

    expect_equal(f$saemParamsToEstimate, c("tka", "tcl", "wt.cl", "tv", "wt.v2", "sexf.cl"))

    expect_equal(setNames(f$saemParHistThetaKeep, f$saemParamsToEstimate),
                 c(tka = 1L, tcl = 1L, wt.cl = 1L, tv = 1L, wt.v2 = 1L, sexf.cl = 0L))

    fit1 <- .nlmixr(mod, dat, "saem", control = saemControlFast)

    expect_equal(fit1$theta["sexf.cl"], c(sexf.cl=1.5))

    # Test literalFix (but it can be fast...)
    currentControl <- saemControlFast
    currentControl$literalFix <- FALSE
    fit1 <-.nlmixr(mod, dat, "saem", control = currentControl)

    expect_equal(fit1$theta["sexf.cl"], c(sexf.cl=1.5))
  })

  test_that("all non time-varying covs", {
    withr::with_options(list(nlmixr2.saemMuRefCov=FALSE), {

      mod <- function() {
        ini({
          tka <- 0.45 ; label("Log Ka")
          tcl <- 1 ; label("Log Cl")
          tv <- 3.45 ; label("Log V")
          eta.ka ~ 0.6
          eta.cl ~ 0.3
          eta.v ~ 0.1
          add.sd <- 0.7
          wt.v2 <- -3
          sexf.cl <- 1.5
          wt.cl <- 3
        })
        model({
          ka <- exp(tka + eta.ka)
          cl <- exp(tcl + eta.cl + wt.cl*logWt70 + sexf.cl*sexf)
          v <- exp(tv + eta.v + wt.v2*logWt70)
          d / dt(depot) <- -ka * depot
          d / dt(center) <- ka * depot - cl / v * center
          ipre <- center / v
          ipre ~ add(add.sd)
        })
      }

      f <- .nlmixr(mod)

      expect_equal(f$saemParamsToEstimateCov, f$saemParamsToEstimate)

      # fit1 <- .nlmixr(mod, dat, "saem", control = saemControlFast))
    })
  })
})
