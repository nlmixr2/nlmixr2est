# saemControl(ebe="mode") for general-likelihood endpoints (dnorm(), t(), ll())
# and censored data: each subject's mode must match FOCEi's posthoc etas at the
# same final parameters, which score the same conditional density.

nmTest({
  .d <- nlmixr2data::theo_sd
  .obs <- .d$EVID == 0
  .dCens <- .d
  .dCens$CENS <- 0
  .dCens$CENS[.obs & .dCens$DV < 2] <- 1
  .dCens$DV[.dCens$CENS == 1] <- 2

  .mod <- function(err, extraIni = "", sdIni = "add.sd <- 0.7") {
    eval(parse(
      text = sprintf(
        "function() {
        ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; %s; %s
          eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1 })
        model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
          cp <- linCmt(); %s })
      }",
        sdIni,
        extraIni,
        err
      )
    ))
  }
  .ctl <- saemControl(print = 0, nBurn = 100, nEm = 100, seed = 42)
  .expectModeIsPosthoc <- function(fit, data) {
    .ph <- .nlmixr(fit, data, est = "focei", control = foceiControl(print = 0, maxOuterIterations = 0L, covMethod = ""))
    .cols <- c("eta.ka", "eta.cl", "eta.v")
    expect_lt(max(abs(as.matrix(fit$eta[, .cols]) - as.matrix(.ph$eta[, .cols]))), 1e-3)
  }

  test_that("the mode of an ll() endpoint is FOCEi's posthoc", {
    skip_on_cran()
    .f <- .nlmixr(
      .mod(
        "sd <- exp(lsd); ll(err) ~ -log(sd) - 0.5 * log(2 * pi) - 0.5 * ((DV - cp) / sd)^2",
        sdIni = "lsd <- log(0.7)"
      ),
      .d,
      est = "saem",
      control = .ctl
    )
    .expectModeIsPosthoc(.f, .d)
  })

  test_that("the normal-model mode uses M2 (uncensored rows with a finite LIMIT)", {
    skip_on_cran()
    # every row carries an M2 term, and some predictions cross LIMIT
    .dM2 <- .d[!.obs | .d$DV > 1, ]
    .dM2$CENS <- 0
    .dM2$LIMIT <- ifelse(.dM2$EVID == 0, 1, NA_real_)
    .f <- .nlmixr(.mod("cp ~ add(add.sd)"), .dM2, est = "saem", control = .ctl)
    .expectModeIsPosthoc(.f, .dM2)
  })

  test_that("censored dnorm() and t() endpoints: saem and its mode use the censoring", {
    skip_on_cran()
    .norm <- .nlmixr(.mod("cp ~ add(add.sd)"), .dCens, est = "saem", control = .ctl)
    .dn <- .nlmixr(.mod("cp ~ add(add.sd) + dnorm()"), .dCens, est = "saem", control = .ctl)
    # the same censored normal likelihood as the closed-form fit
    .th <- c("tka", "tcl", "tv")
    expect_equal(fixef(.dn)[.th], fixef(.norm)[.th], tolerance = 0.02)
    .expectModeIsPosthoc(.dn, .dCens)
    .t <- .nlmixr(.mod("cp ~ add(add.sd) + dt(nu)", "nu <- 5"), .dCens, est = "saem", control = .ctl)
    .expectModeIsPosthoc(.t, .dCens)
  })
})
