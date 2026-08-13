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

  test_that(".updateParFixedRefreshSeFromCov updates parFixedDf and the formatted parFixed (#816)", {
    .uiMod816 <- function() {
      ini({
        tcl <- 1; tfix <- fix(2); add.sd <- 0.7
        eta.cl ~ 0.3
      })
      model({
        cl <- exp(tcl + eta.cl)
        v <- tfix
        d/dt(center) <- -cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .pf <- data.frame(
      Parameter = c("", "", ""),
      Estimate = c(1, 2, 0.7),
      SE = c(0.1, NA_real_, NA_real_),
      "%RSE" = c(10, NA_real_, NA_real_),
      "Back-transformed" = c(exp(1), 2, 0.7),
      "CI Lower" = c(exp(1 - qnorm(0.975) * 0.1), NA_real_, NA_real_),
      "CI Upper" = c(exp(1 + qnorm(0.975) * 0.1), NA_real_, NA_real_),
      check.names = FALSE,
      row.names = c("tcl", "tfix", "add.sd"))
    env <- new.env(parent = emptyenv())
    env$ui <- rxode2::assertRxUi(.uiMod816)
    env$parFixedDf <- .pf
    env$parFixed <- .updateParFixedApplySig(.pf, 3L, 0.95, "tfix")
    class(env$parFixed) <- c("nlmixr2ParFixed", "data.frame")
    .cov <- matrix(c(0.04, 0, 0, 0.0025), 2, 2,
                   dimnames = list(c("tcl", "add.sd"), c("tcl", "add.sd")))

    .updateParFixedRefreshSeFromCov(env, .cov, onlyMissing = TRUE)

    # add.sd filled from sqrt(diag(cov)); tcl kept (onlyMissing = TRUE)
    expect_equal(unname(env$parFixedDf["add.sd", "SE"]), 0.05)
    expect_equal(unname(env$parFixedDf["add.sd", "%RSE"]), 0.05 / 0.7 * 100)
    expect_equal(unname(env$parFixedDf["tcl", "SE"]), 0.1)
    # CI recomputed on the natural scale (Back-transformed == Estimate)
    expect_equal(unname(env$parFixedDf["add.sd", "CI Lower"]),
                 0.7 - qnorm(0.975) * 0.05)
    # formatted table regenerated; FIXED decoration preserved
    expect_equal(as.numeric(env$parFixed["add.sd", "SE"]), 0.05)
    expect_equal(env$parFixed["tfix", "SE"], "FIXED")

    # onlyMissing = FALSE overwrites the structural SE from the cov, and the
    # CI of a log-scale theta is recomputed through its default back-transform
    .updateParFixedRefreshSeFromCov(env, .cov)
    expect_equal(unname(env$parFixedDf["tcl", "SE"]), 0.2)
    expect_equal(as.numeric(env$parFixed["tcl", "SE"]), 0.2)
    expect_equal(unname(env$parFixedDf["tcl", "CI Lower"]),
                 exp(1 - qnorm(0.975) * 0.2))
    expect_equal(unname(env$parFixedDf["tcl", "CI Upper"]),
                 exp(1 + qnorm(0.975) * 0.2))

    # a denormal (uninitialized-memory signature) counts as missing
    env$parFixedDf["add.sd", "SE"] <- 9.39e-323
    .updateParFixedRefreshSeFromCov(env, .cov, onlyMissing = TRUE)
    expect_equal(unname(env$parFixedDf["add.sd", "SE"]), 0.05)

    # A non-default ci lives in the fit's control, not the ui (the ui slot keeps
    # the default), so reading only the ui recomputed the bounds at 95% and
    # labeled the column 95% over an 80% interval.
    env$control <- list(ci = 0.8)
    env$parFixedDf <- .pf
    env$parFixed <- .updateParFixedApplySig(.pf, 3L, 0.8, "tfix")
    class(env$parFixed) <- c("nlmixr2ParFixed", "data.frame")
    .updateParFixedRefreshSeFromCov(env, .cov)
    expect_equal(unname(env$parFixedDf["tcl", "CI Lower"]),
                 exp(1 - qnorm(0.9) * 0.2))
    expect_equal(unname(env$parFixedDf["tcl", "CI Upper"]),
                 exp(1 + qnorm(0.9) * 0.2))
    expect_true("Back-transformed(80%CI)" %in% names(env$parFixed))

    # a control without a usable ci falls back to the ui / the 0.95 default
    env$control <- list(ci = NULL)
    env$parFixedDf <- .pf
    env$parFixed <- .updateParFixedApplySig(.pf, 3L, 0.95, "tfix")
    class(env$parFixed) <- c("nlmixr2ParFixed", "data.frame")
    .updateParFixedRefreshSeFromCov(env, .cov)
    expect_equal(unname(env$parFixedDf["tcl", "CI Lower"]),
                 exp(1 - qnorm(0.975) * 0.2))
    expect_true("Back-transformed(95%CI)" %in% names(env$parFixed))
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
