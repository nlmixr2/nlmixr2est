nmTest({
  test_that("SD is consistent between parFixedDf and parFixed", {
    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- exp(0.45) # Log Ka
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
        ka <- tka + eta.ka
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    fit <- nlmixr(one.cmt, theo_sd, est="focei",
                  control=list(print=0))


    expect_equal(fit$parFixedDf[["BSV(CV% or SD)"]][1],
                 sqrt(fit$omega[1, 1]))

    # parFixed prints the SD to the foceiControl() default `sigdig`
    # significant figures (now 4); parFixedDf keeps full precision.
    expect_equal(as.numeric(fit$parFixed[["BSV(CV% or SD)"]][1]),
                 signif(sqrt(fit$omega[1, 1]), foceiControl()$sigdig))
  })

  ## Literally-fixed thetas are re-inserted into $popDf after the C++ step, so
  ## the default (exp/expit/probitInv) back-transform must be applied in R;
  ## previously they printed on the raw (log/logit) scale.
  test_that("literally-fixed thetas back-transform in the parFixed table", {
    covMod <- function() {
      ini({
        tcl <- 1; beta_cl_wt <- 0.5; tv <- 3.45; add.sd <- 0.7
        eta.cl ~ 0.3
      })
      model({
        cl <- exp(tcl + beta_cl_wt * log(WT / 70) + eta.cl); v <- exp(tv)
        d/dt(center) = -cl / v * center
        cp = center / v
        cp ~ add(add.sd)
      })
    }
    ui <- rxode2::assertRxUi(covMod)
    ## exp() mu-parameter back-transforms
    expect_equal(.updateParFixedBackTransformFixed(ui, "tcl", 1), exp(1))
    ## covariate coefficient and residual error stay on the natural scale
    expect_equal(.updateParFixedBackTransformFixed(ui, "beta_cl_wt", 0.5), 0.5)
    expect_equal(.updateParFixedBackTransformFixed(ui, "add.sd", 0.7), 0.7)

    boundMod <- function() {
      ini({ tf <- 0.5; tv <- 3.45; add.sd <- 0.7; eta.v ~ 0.1 })
      model({
        f <- expit(tf); v <- exp(tv + eta.v)
        d/dt(center) = -center
        cp = f * center / v
        cp ~ add(add.sd)
      })
    }
    uiB <- rxode2::assertRxUi(boundMod)
    ## expit() bounded parameter back-transforms with its bounds
    expect_equal(.updateParFixedBackTransformFixed(uiB, "tf", 0.5),
                 rxode2::expit(0.5))
  })

  ## End-to-end: a fixed mu-referenced log parameter reports exp(value), not the
  ## raw log-scale value, in $parFixedDf.
  test_that("a fixed structural parameter is back-transformed end to end", {
    fixMod <- function() {
      ini({
        tka <- 0.45; tcl <- fix(1); tv <- 3.45; add.sd <- 0.7
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.sd)
      })
    }
    fit <- nlmixr(fixMod, theo_sd, est = "focei",
                  control = foceiControl(print = 0L, maxOuterIterations = 0L,
                                         maxInnerIterations = 0L, calcTables = FALSE))
    expect_equal(fit$parFixedDf["tcl", "Back-transformed"], exp(1),
                 tolerance = 1e-6)
    expect_equal(fit$parFixedDf["tcl", "Estimate"], 1, tolerance = 1e-6)
  })
})
