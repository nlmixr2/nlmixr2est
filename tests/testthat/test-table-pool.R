nmTest({
  # The pooled table (#1356) against the rxSolve() one, across the model and
  # data features the table step handles.
  .poolOne <- function(err = "add") {
    .m <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .m
  }
  .poolFit <- function(model, data, est = "focei", control = NULL) {
    if (is.null(control)) {
      control <- foceiControl(print = 0, maxOuterIterations = 0, calcTables = FALSE)
    }
    suppressMessages(suppressWarnings(nlmixr2(model, data, est, control)))
  }
  .poolExpect <- function(fit, ..., pooled = TRUE) {
    .old <- rxode2::rxWithSeed(
      1356,
      withr::with_options(
        list(nlmixr2.tablePool = FALSE),
        suppressMessages(suppressWarnings(addTable(fit, ...)))
      )
    )
    .i0 <- .tablePoolInfo()
    .new <- rxode2::rxWithSeed(1356, suppressMessages(suppressWarnings(addTable(fit, ...))))
    expect_equal(.tablePoolInfo()$pooledN - .i0$pooledN, as.integer(pooled))
    expect_equal(as.data.frame(.new), as.data.frame(.old), tolerance = 0)
    if (pooled) {
      expect_equal(.odeSwapInfo()$nLoaded, 0L)
    }
    invisible(.new)
  }
  .theo <- nlmixr2data::theo_sd

  test_that("pooled tables: linCmt(), Box-Cox and a lag", {
    .lin <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .poolExpect(.poolFit(.lin, .theo))
    .bc <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
        lambda <- 0.5
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd) + boxCox(lambda)
      })
    }
    .poolExpect(.poolFit(.bc, .theo))
    .lag <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        tlag <- -1
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        alag(depot) <- exp(tlag)
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .poolExpect(.poolFit(.lag, .theo))
  })

  test_that("pooled tables: censoring, EVID=2 rows and a single subject", {
    .d <- .theo
    .d$CENS <- 0
    .d$LIMIT <- NA
    .w <- .d$EVID == 0 & .d$DV < 1
    .d$CENS[.w] <- 1
    .d$DV[.w] <- 1
    .fit <- .poolFit(.poolOne(), .d)
    .poolExpect(.fit)
    # unseeded, the pooled route leaves R's RNG where the rxSolve() route does
    withr::with_seed(3, {
      .a <- suppressMessages(addTable(.fit))
      .ra <- runif(1)
    })
    withr::with_seed(3, {
      .b <- withr::with_options(list(nlmixr2.tablePool = FALSE), suppressMessages(addTable(.fit)))
      .rb <- runif(1)
    })
    expect_equal(as.data.frame(.a), as.data.frame(.b), tolerance = 0)
    expect_equal(.ra, .rb)

    .d2 <- .theo
    .e2 <- .d2[.d2$EVID == 0 & !duplicated(.d2$ID, fromLast = TRUE), ]
    .e2$EVID <- 2
    .e2$TIME <- .e2$TIME + 6
    .e2$DV <- NA
    .d2 <- rbind(.d2, .e2)
    .d2 <- .d2[order(.d2$ID, .d2$TIME), ]
    .poolExpect(.poolFit(.poolOne(), .d2))

    .poolExpect(.poolFit(.poolOne(), .theo[.theo$ID == 1, ]))
  })

  test_that("pooled tables: time-varying covariates, kept columns and npde", {
    .tv <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl) * (WT / 70)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .d <- .theo
    .d$WT <- .d$WT + .d$TIME / 5
    .fit <- .poolFit(.tv, .d)
    .poolExpect(.fit)
    .poolExpect(.fit, table = tableControl(keep = "WT"))
    .poolExpect(.fit, table = tableControl(npde = TRUE))
    .poolExpect(.poolFit(
      .tv,
      .d,
      control = foceiControl(
        print = 0,
        maxOuterIterations = 0,
        calcTables = FALSE,
        rxControl = rxode2::rxControl(covsInterpolation = "linear")
      )
    ))
  })

  test_that("pooled tables: two endpoints, also once the pool is gone", {
    .pkpd <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
        e0 <- 1
        prop.sd <- 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        eff <- e0 + cp
        cp ~ add(add.sd)
        eff ~ prop(prop.sd)
      })
    }
    .d <- .theo
    .d$CMT <- ifelse(.d$EVID == 0, "cp", "depot")
    .e <- .d[.d$EVID == 0, ]
    .e$CMT <- "eff"
    .e$DV <- 1 + .e$DV * 1.1
    .d <- rbind(.d, .e)
    .d <- .d[order(.d$ID, .d$TIME, .d$CMT != "depot"), ]
    .fit <- .poolFit(.pkpd, .d)
    .ref <- .poolExpect(.fit)
    # another fit and a plain rxSolve() in between
    invisible(.poolFit(.poolOne(), .theo))
    invisible(rxode2::rxSolve(.fit, .d))
    .tab <- rxode2::rxWithSeed(1356, suppressMessages(addTable(.fit)))
    expect_equal(as.data.frame(.tab), as.data.frame(.ref), tolerance = 0)
  })

  test_that("pooled tables: saem and posthoc fits", {
    .poolExpect(.poolFit(.poolOne(), .theo, "saem", saemControl(print = 0, nBurn = 10, nEm = 10, calcTables = FALSE)))
    .poolExpect(.poolFit(.poolOne(), .theo, "posthoc", foceiControl(print = 0, calcTables = FALSE)))
  })

  test_that("pooled tables: tolerance-only recovery and a method the pool does not solve", {
    .fit <- .poolFit(.poolOne(), .theo)
    # fails at the fit's atol, recovers once loosened: no method switch, no warning
    withr::with_options(list(nlmixr2.tableFailAtol = .fit$atol[1]), {
      .i0 <- .tablePoolInfo()
      .poolExpect(.fit)
      expect_gt(.tablePoolInfo()$rebuildN, .i0$rebuildN)
    })
    # a fit solved with an ODE method outside the pool's takes the rxSolve() route
    skip_if_not("cvode" %in% names(rxode2::odeMethodToInt(NULL)))
    .fit2 <- .poolFit(
      .poolOne(),
      .theo,
      control = foceiControl(
        print = 0,
        maxOuterIterations = 0,
        calcTables = FALSE,
        rxControl = rxode2::rxControl(method = "cvode")
      )
    )
    .poolExpect(.fit2, pooled = FALSE)
    expect_equal(.tablePoolInfo()$lastDecline, "method")
  })

  # infusions, an addl bolus and a modeled duration, for the dosing rows
  .doseData <- function(nid = 8) {
    .one <- function(id) {
      .d <- data.frame(
        ID = id,
        TIME = c(0, 0, 1, 2, 4, 6, 8, 12, 24, 25, 26, 30, 36, 48),
        EVID = 0,
        AMT = 0,
        RATE = 0,
        ADDL = 0,
        II = 0,
        CMT = 2,
        DV = NA
      )
      .d$EVID[1] <- 1
      .d$AMT[1] <- 100
      .d$RATE[1] <- 50
      .d$EVID[2] <- 1
      .d$AMT[2] <- 80
      .d$CMT[2] <- 1
      .d$ADDL[2] <- 2
      .d$II[2] <- 12
      .d$EVID[9] <- 1
      .d$AMT[9] <- 60
      .d$RATE[9] <- -2
      .d$DV[.d$EVID == 0] <- 2 * exp(rnorm(sum(.d$EVID == 0), 0, 0.2))
      .d
    }
    do.call(rbind, lapply(seq_len(nid), .one))
  }
  .doseModel <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      tdur <- 0.5
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      dur(center) <- exp(tdur)
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  test_that("pooled tables: dosing rows of infusions, addl and a modeled duration", {
    skip_if_not(.tablePoolRxWalk(), "this rxode2 has no output-record walk")
    .d <- rxode2::rxWithSeed(4, .doseData())
    .fit <- .poolFit(.doseModel, .d)
    for (.sn in c(TRUE, FALSE)) {
      .poolExpect(.fit, table = tableControl(addDosing = TRUE, subsetNonmem = .sn))
      .poolExpect(.fit, table = tableControl(addDosing = TRUE, subsetNonmem = .sn, cwres = FALSE))
      .poolExpect(.fit, table = tableControl(addDosing = TRUE, subsetNonmem = .sn, npde = TRUE, nsim = 100))
    }
    # a single subject: no id column from the solve
    .one <- .poolFit(.doseModel, .d[.d$ID == 1, ])
    .poolExpect(.one, table = tableControl(addDosing = TRUE))
    .poolExpect(.one, table = tableControl(addDosing = TRUE, subsetNonmem = FALSE))
  })

  test_that("pooled tables: a mixture, each subject in its own component", {
    .mixMod <- function() {
      ini({
        tka <- 0.45
        tcl1 <- log(c(0, 2.7, 100))
        tcl2 <- log(c(0, 0.1, 120))
        tv <- 3.45
        p1 <- 0.3
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    # a second group with three times the concentrations: the low clearance
    .d2 <- .theo
    .d2$ID <- .d2$ID + 12
    .d2$DV <- .d2$DV * 3
    .fit <- withr::with_seed(
      42,
      .poolFit(
        .mixMod,
        rbind(.theo, .d2),
        control = foceiControl(print = 0, maxOuterIterations = 5, calcTables = FALSE)
      )
    )
    expect_length(unique(.residMixIcov(.fit)$mixest), 2L)
    .tab <- .poolExpect(.fit)
    expect_setequal(unique(.tab$mixest), c(1, 2))
    .poolExpect(.fit, table = tableControl(npde = TRUE, nsim = 100))
    if (.tablePoolRxWalk()) {
      .poolExpect(.fit, table = tableControl(addDosing = TRUE))
    }
  })
})
