# Sensitivities through lag()/diff() of a calculated variable (#1176)
nmTest({
  .lagMod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.cl ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      c0 <- central / v
      cp <- 0.5 * c0 + 0.5 * lag(c0)
      cp ~ add(add.sd)
    })
  }

  # eta also enters directly, a lagged variable of a lagged variable,
  # lag(v, 1), and a residual variance that moves with the lags
  .lagMod2 <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.cl ~ 0.1
      eta.v ~ 0.1
      prop.sd <- 0.1
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      c0 <- central / v
      c1 <- c0 * exp(eta.v)
      cp <- 0.3 * c0 + 0.4 * lag(c0, 1) + 0.3 * lag(c1)
      cp ~ prop(prop.sd)
    })
  }

  # a lagged variable inside an ODE, and diff()
  .lagMod3 <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      tke <- -1
      eta.cl ~ 0.1
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      ke <- exp(tke)
      c0 <- central / v
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      d/dt(eff) <- ke * (c0 - eff)
      cp <- eff + 0.5 * diff(c0)
      cp ~ add(add.sd)
    })
  }

  .lagDat <- nlmixr2data::theo_sd[nlmixr2data::theo_sd$ID == 1, ]

  # max |analytic - central difference| of rx__sens_<v>_BY_<par>___ per
  # parameter, from a solve of the model text `txt`
  .lagFd <- function(mod, txt, par, idx, at = 0.2) {
    .ui <- rxode2::rxode2(mod)
    .m <- rxode2::rxode2(txt)
    .th <- .ui$iniDf$est[!is.na(.ui$iniDf$ntheta)]
    .neta <- max(.ui$iniDf$neta1, na.rm = TRUE)
    .p <- c(
      setNames(.th, paste0("THETA[", seq_along(.th), "]")),
      setNames(rep(at, .neta), paste0("ETA[", seq_len(.neta), "]"))
    )
    .sol <- function(p) {
      as.data.frame(rxode2::rxSolve(.m, p, .lagDat, atol = 1e-12, rtol = 1e-12, addDosing = FALSE))
    }
    .s0 <- .sol(.p)
    .h <- 1e-5
    .ret <- NULL
    for (k in idx) {
      .nm <- paste0(par, "[", k, "]")
      .pp <- .p
      .pp[.nm] <- .p[.nm] + .h
      .pm <- .p
      .pm[.nm] <- .p[.nm] - .h
      .sp <- .sol(.pp)
      .sm <- .sol(.pm)
      for (v in c("rx_pred_", "rx_r_")) {
        .fd <- (.sp[[v]] - .sm[[v]]) / (2 * .h)
        .an <- .s0[[paste0("rx__sens_", v, "_BY_", par, "_", k, "___")]]
        .ret <- rbind(.ret, data.frame(k = k, v = v, err = max(abs(.fd - .an)), fd = max(abs(.fd))))
      }
    }
    .ret
  }

  test_that("only lag() of a calculated variable needs the chain", {
    expect_true(.foceiUsesLagVar(rxode2::rxode2(.lagMod)))
    expect_true(.foceiUsesLagVar(rxode2::rxode2(.lagMod3)))
    .modWt <- function() {
      ini({
        tcl <- 1
        tv <- 3.45
        eta.cl ~ 0.1
        add.sd <- 0.7
      })
      model({
        cl <- exp(tcl + eta.cl) * lag(WT) / 70
        v <- exp(tv)
        lag(central) <- 0.1
        d/dt(central) <- -cl / v * central
        cp <- central / v
        cp ~ add(add.sd)
      })
    }
    expect_false(.foceiUsesLagVar(rxode2::rxode2(.modWt)))
  })

  test_that("inner eta sensitivities chain through lag()", {
    for (.mod in list(.lagMod, .lagMod2, .lagMod3)) {
      .s <- rxode2::rxode2(.mod)$foceiEnv
      expect_true(length(.s$..lagSens) > 0L)
      .fd <- .lagFd(.mod, .s$..inner, "ETA", seq_len(.s$..maxEta))
      expect_true(all(.fd$fd[.fd$v == "rx_pred_"] > 0.1))
      expect_true(all(.fd$err < 1e-6 * pmax(1, .fd$fd)))
    }
  })

  test_that("impmap theta sensitivities chain through lag()", {
    .ts <- rxode2::rxode2(.lagMod2)$impmapThetaSens
    .fd <- .lagFd(.lagMod2, .ts$thetaSens, "THETA", .ts$thetaSensIdx)
    expect_true(any(.fd$fd[.fd$v == "rx_pred_"] > 0.1))
    expect_true(all(.fd$err < 1e-6 * pmax(1, .fd$fd)))
    # the combined eta+theta inner model (#958)
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.lagMod2))
    assign("control", impmapControl(), envir = .ui)
    .s <- .ui$foceiEnv
    expect_true(length(.s$..combThetaIdx) > 0L)
    .fd <- .lagFd(.lagMod2, .s$..inner, "THETA", .s$..combThetaIdx)
    expect_true(all(.fd$err < 1e-6 * pmax(1, .fd$fd)))
  })

  test_that("focei fits a model whose eta reaches the prediction only through lag()", {
    .fit <- nlmixr2(.lagMod, nlmixr2data::theo_sd, est = "focei", control = foceiControl(print = 0))
    expect_true(is.finite(.fit$objf))
    expect_true(sd(.fit$eta$eta.cl) > 0.01)
    expect_message(
      .fast <- nlmixr2(.lagMod, nlmixr2data::theo_sd, est = "focei", control = foceiControl(print = 0, fast = TRUE)),
      "lag\\(\\) of a calculated variable"
    )
    expect_equal(.fast$objf, .fit$objf, tolerance = 1e-6)
  })
})
