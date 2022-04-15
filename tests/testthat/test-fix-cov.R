test_that("test mu reference covariate in saem", {

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

  f1 <- nlmixr(mod)

  f2 <- nlmixr(mod2)
  expect_equal(f1$saemParamsToEstimate, f2$saemParamsToEstimate)
  expect_equal(f1$saemInitTheta, f2$saemInitTheta)


  nid <- 12

  dat <- rxode2::rxWithSeed(102478, {

    cov <- data.frame(id=seq(nid), logWt70=log(rnorm(nid, 70, 10) / 70), sexf=round(runif(nid)))

    ev <- et() %>%
      et(amt=320) %>%
      et(c(0.25, 0.5, 1, 2, 4, 6, 8, 10, 12, 24)) %>%
      et(id=1:nid) %>%
      merge(cov)

    rxSolve(mod, ev, addDosing=TRUE, returnType="data.frame") %>%
      dplyr::select(id, evid, cmt, amt, time, sim, logWt70, sexf) %>%
      dplyr::rename(dv=sim)

  })


  fit1 <- nlmixr(mod, dat, "saem")

  fit2 <- nlmixr(mod2, dat, "saem")

  for (n in names(fit1$theta)) {
    expect_equal(fit1$theta[n], fit2$theta[n])
  }

  expect_false(all(fit1$theta == fit2$theta))

})
