nmTest({
  test_that("FOCEi HdEta emits linCmtB(which1=-3) for a modeled alag() on a linCmt() compartment (#920)", {
    one.cmt.lag <- function() {
      ini({
        tka <- log(1.15)
        tcl <- log(0.135)
        tv <- log(8)
        tlag <- log(0.2)
        eta.ka ~ 0.5
        eta.cl ~ 0.1
        eta.v ~ 0.1
        eta.lag ~ 0.5
        prop.err <- 0.15
        add.err <- 0.6
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        tlagi <- exp(tlag + eta.lag)
        alag(depot) <- tlagi
        linCmt() ~ prop(prop.err) + add(add.err)
      })
    }

    ui <- one.cmt.lag()
    s <- rxUiGet.foceiEnv(list(ui))

    # eta.lag is ETA_4_ -- before #920 this line rendered as the literal "0"
    .lagLine <- s$..HdEta[grepl("BY_ETA_4___", s$..HdEta, fixed = TRUE)]
    expect_true(any(grepl("linCmtB(", .lagLine, fixed = TRUE)))
    expect_true(any(grepl("-3", .lagLine, fixed = TRUE)))
    expect_false(any(grepl("^rx__sens_rx_pred__BY_ETA_4___=0$", .lagLine)))

    # rx_r_ (proportional error) embeds rx_pred_'s linCmtB() call directly, so
    # its own d(rx_r_)/d(eta.lag) needs the same correction chained through it
    .lagLineR <- s$..REta[grepl("BY_ETA_4___", s$..REta, fixed = TRUE)]
    expect_true(any(grepl("linCmtB(", .lagLineR, fixed = TRUE)))
    expect_false(any(grepl("^rx__sens_rx_r__BY_ETA_4___=0$", .lagLineR)))
  })

  test_that("FOCEi HdEta/REta emit the pred/F correction for a modeled f() on a linCmt() compartment (#920)", {
    one.cmt.f <- function() {
      ini({
        tka <- log(1.15)
        tcl <- log(0.135)
        tv <- log(8)
        tf <- log(0.8)
        eta.ka ~ 0.5
        eta.cl ~ 0.1
        eta.v ~ 0.1
        eta.f ~ 0.3
        prop.err <- 0.15
        add.err <- 0.6
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        fi <- exp(tf + eta.f)
        f(depot) <- fi
        linCmt() ~ prop(prop.err) + add(add.err)
      })
    }

    ui <- one.cmt.f()
    s <- rxUiGet.foceiEnv(list(ui))

    .fLine <- s$..HdEta[grepl("BY_ETA_4___", s$..HdEta, fixed = TRUE)]
    expect_false(any(grepl("^rx__sens_rx_pred__BY_ETA_4___=0$", .fLine)))

    # rx_r_ needs the same chained correction as the alag() case above
    .fLineR <- s$..REta[grepl("BY_ETA_4___", s$..REta, fixed = TRUE)]
    expect_false(any(grepl("^rx__sens_rx_r__BY_ETA_4___=0$", .fLineR)))
  })

  test_that("pred/F correction matches a central finite difference (#920)", {
    one.cmt.f <- function() {
      ini({
        tka <- log(1.15)
        tcl <- log(0.135)
        tv <- log(8)
        tf <- log(0.8)
        eta.ka ~ 0.5
        eta.cl ~ 0.1
        eta.v ~ 0.1
        eta.f ~ 0.3
        prop.err <- 0.15
        add.err <- 0.6
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        fi <- exp(tf + eta.f)
        f(depot) <- fi
        linCmt() ~ prop(prop.err) + add(add.err)
      })
    }

    ui <- one.cmt.f()
    s <- rxUiGet.foceiEnv(list(ui))
    nlmixr2est:::.rxFinalizeInner(s, FALSE, FALSE, 0L)
    mod <- suppressMessages(rxode2::rxode2(s$..inner))

    ev <- rxode2::et(amt = 100, cmt = "depot") %>% rxode2::et(seq(0.1, 24, by = 0.5))
    THETA <- c(log(1.15), log(0.135), log(8), log(0.8), 0.15, 0.6)
    ETA0 <- c(0.1, -0.05, 0.02, 0.2)

    solveAt <- function(eta4) {
      ETA <- ETA0
      ETA[4] <- eta4
      pars <- c(stats::setNames(THETA, paste0("THETA[", 1:6, "]")),
                stats::setNames(ETA, paste0("ETA[", 1:4, "]")))
      suppressMessages(rxode2::rxSolve(mod, pars, ev, returnType = "data.frame"))
    }

    h <- 1e-4
    sPlus <- solveAt(ETA0[4] + h)
    sMinus <- solveAt(ETA0[4] - h)
    s0 <- solveAt(ETA0[4])

    fdPred <- (sPlus$rx_pred_ - sMinus$rx_pred_) / (2 * h)
    anaPred <- s0$rx__sens_rx_pred__BY_ETA_4___
    expect_false(anyNA(anaPred))
    expect_equal(anaPred, fdPred, tolerance = 1e-5)

    fdR <- (sPlus$rx_r_ - sMinus$rx_r_) / (2 * h)
    anaR <- s0$rx__sens_rx_r__BY_ETA_4___
    expect_equal(anaR, fdR, tolerance = 1e-5)
  })

  test_that("more than one f()-scaled linCmt() compartment leaves the gradient uncorrected (rxode2/rxode2#1237)", {
    two.f <- function() {
      ini({
        tka <- log(1.15)
        tcl <- log(0.135)
        tv <- log(8)
        tq <- log(1)
        tvp <- log(10)
        tf1 <- log(0.7)
        tf2 <- log(0.9)
        eta.ka ~ 0.5
        eta.cl ~ 0.1
        eta.v ~ 0.1
        eta.f1 ~ 0.2
        eta.f2 ~ 0.2
        prop.err <- 0.15
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        q <- exp(tq)
        vp <- exp(tvp)
        f1 <- exp(tf1 + eta.f1)
        f2 <- exp(tf2 + eta.f2)
        f(depot) <- f1
        f(central) <- f2
        linCmt() ~ prop(prop.err)
      })
    }

    ui <- two.f()
    s <- suppressWarnings(rxUiGet.foceiHdEta(list(ui)))

    expect_true(any(grepl("^rx__sens_rx_pred__BY_ETA_4___=0$", s$..HdEta)))
    expect_true(any(grepl("^rx__sens_rx_pred__BY_ETA_5___=0$", s$..HdEta)))
  })

  test_that("foceiControl(eventSens='fd') opts out of the linCmt() f() correction (rxode2/rxode2#1236)", {
    one.cmt.f <- function() {
      ini({
        tcl <- log(0.135)
        tv <- log(8)
        tf <- log(0.8)
        eta.cl ~ 0.1
        eta.v ~ 0.1
        eta.f ~ 0.3
        prop.err <- 0.15
      })
      model({
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        fi <- exp(tf + eta.f)
        f(central) <- fi
        linCmt() ~ prop(prop.err)
      })
    }

    ui <- rxode2::rxUiDecompress(one.cmt.f())
    rxode2::rxAssignControlValue(ui, "eventSens", "fd")
    s <- suppressWarnings(rxUiGet.foceiHdEta(list(ui)))
    expect_true(any(grepl("^rx__sens_rx_pred__BY_ETA_3___=0$", s$..HdEta)))
  })

  test_that("linCmtB(which1=-3) alag correction matches a central finite difference (#920)", {
    one.cmt.lag <- function() {
      ini({
        tka <- log(1.15)
        tcl <- log(0.135)
        tv <- log(8)
        tlag <- log(0.2)
        eta.ka ~ 0.5
        eta.cl ~ 0.1
        eta.v ~ 0.1
        eta.lag ~ 0.5
        prop.err <- 0.15
        add.err <- 0.6
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        tlagi <- exp(tlag + eta.lag)
        alag(depot) <- tlagi
        linCmt() ~ prop(prop.err) + add(add.err)
      })
    }

    ui <- one.cmt.lag()
    s <- rxUiGet.foceiEnv(list(ui))
    nlmixr2est:::.rxFinalizeInner(s, FALSE, FALSE, 0L)
    mod <- suppressMessages(rxode2::rxode2(s$..inner))

    ev <- rxode2::et(amt = 100, cmt = "depot") %>% rxode2::et(seq(0.1, 24, by = 0.5))
    THETA <- c(log(1.15), log(0.135), log(8), log(0.2), 0.15, 0.6)
    ETA0 <- c(0.1, -0.05, 0.02, 0.3)

    solveAt <- function(eta4) {
      ETA <- ETA0
      ETA[4] <- eta4
      pars <- c(stats::setNames(THETA, paste0("THETA[", 1:6, "]")),
                stats::setNames(ETA, paste0("ETA[", 1:4, "]")))
      suppressMessages(rxode2::rxSolve(mod, pars, ev, returnType = "data.frame"))
    }

    h <- 1e-4
    sPlus <- solveAt(ETA0[4] + h)
    sMinus <- solveAt(ETA0[4] - h)
    s0 <- solveAt(ETA0[4])

    fdPred <- (sPlus$rx_pred_ - sMinus$rx_pred_) / (2 * h)
    anaPred <- s0$rx__sens_rx_pred__BY_ETA_4___
    expect_false(anyNA(anaPred))
    expect_equal(anaPred, fdPred, tolerance = 1e-5)

    fdR <- (sPlus$rx_r_ - sMinus$rx_r_) / (2 * h)
    anaR <- s0$rx__sens_rx_r__BY_ETA_4___
    expect_equal(anaR, fdR, tolerance = 1e-5)
  })

  test_that("more than one lagged linCmt() compartment leaves the gradient uncorrected (rxode2/rxode2#1237)", {
    two.lag <- function() {
      ini({
        tka <- log(1.15)
        tcl <- log(0.135)
        tv <- log(8)
        tq <- log(1)
        tvp <- log(10)
        tlag1 <- log(0.2)
        tlag2 <- log(0.4)
        eta.ka ~ 0.5
        eta.cl ~ 0.1
        eta.v ~ 0.1
        eta.lag1 ~ 0.3
        eta.lag2 ~ 0.3
        prop.err <- 0.15
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        q <- exp(tq)
        vp <- exp(tvp)
        lag1 <- exp(tlag1 + eta.lag1)
        lag2 <- exp(tlag2 + eta.lag2)
        alag(depot) <- lag1
        alag(central) <- lag2
        linCmt() ~ prop(prop.err)
      })
    }

    ui <- two.lag()
    s <- suppressWarnings(rxUiGet.foceiHdEta(list(ui)))

    expect_true(any(grepl("^rx__sens_rx_pred__BY_ETA_4___=0$", s$..HdEta)))
    expect_true(any(grepl("^rx__sens_rx_pred__BY_ETA_5___=0$", s$..HdEta)))
  })

  test_that("foceiControl(eventSens='fd') opts out of the linCmt() alag correction (rxode2/rxode2#1236)", {
    one.cmt.lag <- function() {
      ini({
        tcl <- log(0.135)
        tv <- log(8)
        tlag <- log(0.3)
        eta.cl ~ 0.1
        eta.v ~ 0.1
        eta.lag ~ 0.4
        prop.err <- 0.15
      })
      model({
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        tlagi <- exp(tlag + eta.lag)
        alag(central) <- tlagi
        linCmt() ~ prop(prop.err)
      })
    }

    ui <- rxode2::rxUiDecompress(one.cmt.lag())
    rxode2::rxAssignControlValue(ui, "eventSens", "fd")
    s <- suppressWarnings(rxUiGet.foceiHdEta(list(ui)))
    expect_true(any(grepl("^rx__sens_rx_pred__BY_ETA_3___=0$", s$..HdEta)))
  })

  test_that("a modeled alag() on a linCmt() compartment fits with FOCEi (#920)", {
    one.cmt.lag <- function() {
      ini({
        tka <- log(1.57)
        tcl <- log(2.72)
        tv <- log(31.5)
        tlag <- log(0.5)
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        eta.lag ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        lag <- exp(tlag + eta.lag)
        alag(depot) <- lag
        linCmt() ~ add(add.sd)
      })
    }

    fit <- .nlmixr(one.cmt.lag, nlmixr2data::theo_sd, est = "focei",
                   control = foceiControl(maxOuterIterations = 5, maxInnerIterations = 50,
                                          covMethod = "", calcTables = FALSE, print = 0))
    expect_true(is.finite(fit$objDf$OBJF))
  })
})
