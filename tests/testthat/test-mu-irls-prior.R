nmTest({
  test_that("irls mu M-step honours a theta prior (row-scaled weights, #1355)", {
    skip_on_cran()
    skip_if_not(
      exists("rxPriorBuildSpec", envir = asNamespace("rxode2"), inherits = FALSE),
      "rxode2 without the shared prior kernel (nlmixr2/rxode2#1270)"
    )
    theoSd2 <- nlmixr2data::theo_sd
    theoSd2$logWT <- log(theoSd2$WT / 70)
    mod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        allo.cl <- 0.75
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + allo.cl * logWT)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    modPrior <- mod |> ini(prior(tcl) ~ dnorm(0.5, 0.01))
    fit0 <- .nlmixr(mod, theoSd2, "ifocei", ifoceiControl(print = 0))
    fit1 <- .nlmixr(modPrior, theoSd2, "ifocei", ifoceiControl(print = 0))
    expect_true(all(is.finite(fit0$theta)))
    expect_true(all(is.finite(fit1$theta)))
    # the tight prior pulls tcl from its data estimate (~1) to ~0.5
    expect_gt(fit0$theta["tcl"], 0.8)
    expect_equal(unname(fit1$theta["tcl"]), 0.5, tolerance = 0.05)
  })
})
