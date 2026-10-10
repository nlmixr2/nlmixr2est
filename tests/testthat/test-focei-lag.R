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

  # an ODE using a lagged variable defined through another one
  .lagMod4 <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.cl ~ 0.1
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      c0 <- central / v
      c1 <- c0 * exp(eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      d/dt(eff) <- c1 - eff
      cp <- eff + lag(c0) + lag(c1)
      cp ~ add(add.sd)
    })
  }

  .lagDat <- nlmixr2data::theo_sd[nlmixr2data::theo_sd$ID == 1, ]

  # max |analytic - central difference| of rx__sens_<v>_BY_<par>___ per
  # parameter, from a solve of the model text `txt`
  .lagFd <- function(mod, txt, par, idx, at = 0.2, dat = .lagDat) {
    .ui <- rxode2::rxode2(mod)
    .m <- rxode2::rxode2(txt)
    .th <- .ui$iniDf$est[!is.na(.ui$iniDf$ntheta)]
    .neta <- max(.ui$iniDf$neta1, na.rm = TRUE)
    .p <- c(
      setNames(.th, paste0("THETA[", seq_along(.th), "]")),
      setNames(rep(at, .neta), paste0("ETA[", seq_len(.neta), "]"))
    )
    .sol <- function(p) {
      as.data.frame(rxode2::rxSolve(.m, p, dat, atol = 1e-12, rtol = 1e-12, addDosing = FALSE))
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
    for (.mod in list(.lagMod, .lagMod2, .lagMod3, .lagMod4)) {
      .s <- rxode2::rxode2(.mod)$foceiEnv
      expect_true(length(.s$..lagSens) > 0L)
      .fd <- .lagFd(.mod, .s$..inner, "ETA", seq_len(.s$..maxEta))
      expect_true(all(.fd$fd[.fd$v == "rx_pred_"] > 0.1))
      expect_true(all(.fd$err < 1e-6 * pmax(1, .fd$fd)))
    }
  })

  test_that("an ODE can use a lagged variable defined by if/else, not its lag()", {
    .ui <- rxode2::rxode2(.lagMod3)
    .ifOde <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        tke <- -1
        eta.cl ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        ke <- exp(tke)
        if (WT > 70) {
          c0 <- central / v
        } else {
          c0 <- 2 * central / v
        }
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        d/dt(eff) <- ke * (c0 - eff)
        cp <- eff + lag(c0)
        cp ~ add(add.sd)
      })
    }
    .s <- rxode2::rxode2(.ifOde)$foceiEnv
    .fd <- .lagFd(.ifOde, .s$..inner, "ETA", seq_len(.s$..maxEta))
    expect_true(all(.fd$err < 1e-6 * pmax(1, .fd$fd)))
    .histOde <- rxode2::model(.ui, d/dt(eff) <- ke * (lag(c0) - eff))
    expect_error(rxode2::rxode2(.histOde)$foceiEnv, "inside an ODE is not supported")
  })

  test_that("a lagged variable read before it is reassigned keeps its sensitivity", {
    .reMod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.cl ~ 0.1
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        c0 <- central / v
        y <- 2 * c0
        c0 <- c0 * exp(eta.cl)
        cp <- y + lag(c0)
        cp ~ add(add.sd)
      })
    }
    .s <- rxode2::rxode2(.reMod)$foceiEnv
    # rxode2 builds before rxode2#1435 read the final c0 in y
    skip_if_not(grepl("rx_lagv1_c0=c0", .s$..inner, fixed = TRUE))
    .fd <- .lagFd(.reMod, .s$..inner, "ETA", seq_len(.s$..maxEta))
    expect_true(all(.fd$err < 1e-6 * pmax(1, .fd$fd)))
  })

  test_that("an ODE between two assignments of a lagged variable reads the snapshot", {
    # rxode2 builds before rxode2#1445 bind the final c0 in the ODE
    skip_if_not(grepl(
      "rx_lagv1_c0",
      rxode2::rxS("c0=central/10\nd/dt(central)=-c0\nc0=3*c0\ncp=lag(c0)")$..ddt,
      fixed = TRUE
    ))
    .odeMod <- function() {
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
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        c0 <- central / v
        d/dt(eff) <- ke * (c0 - eff)
        c0 <- c0 * exp(eta.cl)
        cp <- eff + lag(c0)
        cp ~ add(add.sd)
      })
    }
    .ui <- rxode2::rxode2(.odeMod)
    .s <- .ui$foceiEnv
    .fd <- .lagFd(.odeMod, .s$..inner, "ETA", seq_len(.s$..maxEta))
    expect_true(all(.fd$err < 1e-6 * pmax(1, .fd$fd)))
    # the inner model predicts what rxSolve() gives for the original model
    .th <- .ui$iniDf[!is.na(.ui$iniDf$ntheta), ]
    .inner <- rxode2::rxSolve(
      rxode2::rxode2(.s$..inner),
      c(setNames(.th$est, paste0("THETA[", seq_along(.th$est), "]")), "ETA[1]" = 0.2, "ETA[2]" = 0.2),
      .lagDat,
      atol = 1e-12,
      rtol = 1e-12,
      addDosing = FALSE
    )
    .orig <- rxode2::rxSolve(
      .ui$simulationModel,
      c(setNames(.th$est, .th$name), eta.cl = 0.2, eta.v = 0.2),
      .lagDat,
      atol = 1e-12,
      rtol = 1e-12,
      addDosing = FALSE
    )
    expect_equal(.inner$rx_pred_, .orig$cp, tolerance = 1e-8)
  })

  test_that("the linCmt() sensitivity carry keeps the lag() terms", {
    .carryMod <- function() {
      ini({
        tka <- 0.45
        tcl <- log(2)
        tv <- log(20)
        eta.cl ~ 0.1
        add.sd <- 0.5
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
        v <- exp(tv)
        c0 <- exp(eta.cl) * wt
        cp <- linCmt() + 0.01 * lag(c0)
        cp ~ add(add.sd)
      })
    }
    .s <- rxode2::rxode2(.carryMod)$foceiEnv
    expect_false(is.null(.s$..linCmtCarryPairs))
    # a time-varying covariate on a linCmt() parameter
    .dat <- within(.lagDat, wt <- WT * (1 + 0.02 * TIME))
    .fd <- .lagFd(.carryMod, .s$..inner, "ETA", 1L, dat = .dat)
    expect_true(all(.fd$err < 1e-6 * pmax(1, .fd$fd)))
    # a residual variance reaching the eta through the lagged variable too
    .propMod <- rxode2::model(rxode2::rxode2(.carryMod), cp ~ prop(add.sd * c0))
    .s <- rxode2::rxode2(.propMod)$foceiEnv
    expect_false(is.null(.s$..linCmtCarryPairs))
    .fd <- .lagFd(.propMod, .s$..inner, "ETA", 1L, dat = .dat)
    expect_true(.fd$fd[.fd$v == "rx_r_"] > 1e-3)
    expect_true(all(.fd$err < 1e-6 * pmax(1, .fd$fd)))
  })

  test_that("the prediction can be a lag() alone", {
    .lagOnly <- rxode2::rxode2(.lagMod)
    .lagOnly <- rxode2::model(.lagOnly, cp <- lag(c0))
    .s <- rxode2::rxode2(.lagOnly)$foceiEnv
    .fd <- .lagFd(.lagOnly, .s$..inner, "ETA", 1L)
    expect_true(.fd$fd[.fd$v == "rx_pred_"] > 0.1)
    expect_true(all(.fd$err < 1e-6 * pmax(1, .fd$fd)))
  })

  test_that("the AR(1) correction chains through lag()", {
    .arMod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.cl ~ 0.1
        add.sd <- 0.7
        ar1.cor <- 0.5
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        c0 <- central / v
        cp <- 0.5 * c0 + 0.5 * lag(c0)
        cp ~ add(add.sd) + ar(ar1.cor)
      })
    }
    # the norm form rxUiGet.focei() builds
    nlmixr2global$rxArNorm <- TRUE
    on.exit(nlmixr2global$rxArNorm <- FALSE)
    .s <- rxode2::rxode2(.arMod)$foceiEnv
    nlmixr2global$rxArNorm <- FALSE
    expect_true(length(.s$..arEtaSens) > 0L)
    # the AR(1) correction reads rx_time_pk, defined in the model prologue (#1167)
    .fd <- .lagFd(.arMod, .addMtimeLines(.s$..inner, .s), "ETA", 1L)
    expect_true(all(.fd$err < 1e-6 * pmax(1, .fd$fd)))
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
    # no 2nd-order sensitivities through a lagged variable
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.lagMod))
    assign("control", list(innerHessian = "conditional"), envir = .ui)
    expect_error(.ui$foceiEnv, "lag\\(\\) of a variable")
    expect_message(
      nlmixr2(.lagMod, nlmixr2data::theo_sd, est = "focei", control = foceiControl(print = 0, covMethod = "analytic")),
      "lag\\(\\) of a calculated variable"
    )
  })

  test_that("rxode2's rx_lagv snapshot lines are taken with the lagged definitions", {
    .lhs <- c(
      "c0=central/v",
      "rx_lagv1_c0=c0",
      "c0=2*rx_lagv1_c0",
      "cp=eff+lag(c0)",
      "rx_lagv1_cp=cp",
      "c00=1",
      "rx_lagvx_c0=1"
    )
    expect_identical(.foceiIsLagDef(.lhs, "c0"), c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
    .s <- new.env()
    .s$..laggedVars <- c("c0", "rx_ar1")
    .s$..lhs <- c(.lhs, "rx_ar1=1", "rx_lagv1_rx_ar1=rx_ar1")
    expect_identical(.foceiLagDefs(.s), .lhs[1:3])
  })
})
