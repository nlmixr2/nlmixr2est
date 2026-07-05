nmTest({

  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  .mkNlmUi <- function(ctl) {
    .ui <- rxode2::rxUiDecompress(rxode2::as.rxUi(one.cmt))
    assign("control", ctl, envir = .ui)
    class(.ui) <- c("nlm", class(.ui))
    .ui
  }

  .sensStates <- function(m) {
    grep("^rx__sens_.*BY_THETA", rxode2::rxState(m), value = TRUE)
  }

  test_that(".nlmAdjointResolve maps base methods to their s-variants", {
    r <- function(m, sm = "adjoint") {
      nlmixr2est:::.nlmAdjointResolve(
        .mkNlmUi(nlmControl(sensMethod = sm,
                            rxControl = rxode2::rxControl(method = m))))
    }
    expect_false(r("liblsoda", "forward")$useAdjoint)

    expect_equal(r("liblsoda")$sMethodName, "liblsodaadj")
    expect_false(r("liblsoda")$stiff)

    expect_equal(r("rk4")$sMethodName, "rk4s")
    expect_false(r("rk4")$stiff)

    expect_equal(r("ros4")$sMethodName, "ros4s")
    expect_true(r("ros4")$stiff)

    expect_equal(r("dop853")$sMethodName, "dop853s")
  })

  test_that("adjoint thetaGrad emits the same sensitivity columns as forward", {
    fwd <- .mkNlmUi(nlmControl(sensMethod = "forward"))$nlmSensModel$thetaGrad
    fwdStates <- .sensStates(fwd)
    expect_true(length(fwdStates) > 0)
    for (m in c("rk4", "ros4")) {
      adj <- .mkNlmUi(nlmControl(sensMethod = "adjoint",
                                 rxControl = rxode2::rxControl(method = m)))$nlmSensModel$thetaGrad
      expect_setequal(.sensStates(adj), fwdStates)
      ## the adjoint model must carry the transpose-Jacobian sweep lhs
      expect_true(any(grepl("^rx__adjFX", rxode2::rxLhs(adj))))
    }
  })

  test_that("nlm adjoint fit recovers the forward estimates and objective", {
    .d <- nlmixr2data::theo_sd
    .rf <- function(sm, method = NULL) {
      .ctl <- if (is.null(method)) {
        nlmControl(sensMethod = sm, solveType = "grad", print = 0)
      } else {
        nlmControl(sensMethod = sm, solveType = "grad", print = 0,
                   rxControl = rxode2::rxControl(method = method))
      }
      .nlmixr(one.cmt, .d, est = "nlm", .ctl)
    }
    .est <- function(f) setNames(f$parFixedDf$Estimate, rownames(f$parFixedDf))

    fwd <- .rf("forward")
    adjExplicit <- .rf("adjoint", "rk4")   # rk4s
    adjStiff <- .rf("adjoint", "ros4")     # ros4s

    expect_equal(adjExplicit$objf, fwd$objf, tolerance = 1e-4)
    expect_equal(adjStiff$objf, fwd$objf, tolerance = 1e-4)

    expect_equal(.est(adjExplicit), .est(fwd), tolerance = 1e-4)
    expect_equal(.est(adjStiff), .est(fwd), tolerance = 1e-4)
  })

  test_that("nlminb (gradient family) shares the adjoint wiring", {
    .d <- nlmixr2data::theo_sd
    fwd <- .nlmixr(one.cmt, .d, est = "nlminb",
                   nlminbControl(sensMethod = "forward", solveType = "grad", print = 0))
    adj <- .nlmixr(one.cmt, .d, est = "nlminb",
                   nlminbControl(sensMethod = "adjoint", solveType = "grad", print = 0,
                                 rxControl = rxode2::rxControl(method = "ros4")))
    expect_equal(adj$objf, fwd$objf, tolerance = 1e-4)
    expect_equal(setNames(adj$parFixedDf$Estimate, rownames(adj$parFixedDf)),
                 setNames(fwd$parFixedDf$Estimate, rownames(fwd$parFixedDf)),
                 tolerance = 1e-4)
  })

  test_that("focei adjoint inner sensitivities match the forward fit", {
    .d <- nlmixr2data::theo_sd
    ## 3 etas, 2 ODE states; the inner sens is requested as adjoint below
    mod <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.1; eta.cl ~ 0.1; eta.v ~ 0.1
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
    ## the adjoint inner model keeps FOCEi's lhs offsets (pred at lhs[0]) and
    ## appends the rx__adj* sweep lhs after them
    ui <- rxode2::rxUiDecompress(rxode2::as.rxUi(mod))
    assign("control", foceiControl(sensMethod = "adjoint"), envir = ui)
    .lhs <- rxode2::rxModelVars(ui$focei$inner)$lhs
    expect_identical(.lhs[1], "rx_pred_")
    expect_true(any(grepl("^rx__adjFX", .lhs)))
    expect_gt(which(.lhs == "rx__adjFX_0_0__"), 8L)

    ## Inner (eta) equivalence at fixed thetas is the tight check of the adjoint
    ## sensitivities.  The full outer fit is only equal to forward within
    ## optimization tolerance: FOCEi's numeric outer gradient amplifies the tiny
    ## (~1e-5) inner-objective round-off difference between the forward and
    ## adjoint solvers, so the converged fits agree to ~1e-2, not exactly.
    .rf <- function(sm, mo) {
      .nlmixr(mod, .d, est = "focei",
              foceiControl(sensMethod = sm, print = 0, maxOuterIterations = mo,
                           covMethod = ""))
    }
    fwd0 <- .rf("forward", 0L)
    adj0 <- .rf("adjoint", 0L)
    expect_equal(adj0$objf, fwd0$objf, tolerance = 1e-3)

    fwd <- .rf("forward", 30L)
    adj <- .rf("adjoint", 30L)
    expect_equal(adj$objf, fwd$objf, tolerance = 0.05)
    expect_equal(setNames(adj$parFixedDf$Estimate, rownames(adj$parFixedDf)),
                 setNames(fwd$parFixedDf$Estimate, rownames(fwd$parFixedDf)),
                 tolerance = 0.05)
  })

  test_that("focei adjoint handles modeled dosing modifiers (lag)", {
    ## 4 etas > 2 states -> adjoint; a modeled lag() exercises the adjoint
    ## dose-jump path (rx__adjDlag).  The full outer fit of this warfarin model
    ## is numerically unstable across sensitivity methods (forward "jump" and
    ## "fd" converge to different optima too), so equivalence is checked at fixed
    ## thetas (inner eta optimization only), where the adjoint dose-jump matches
    ## the analytic-jump forward sensitivities.
    lag.mod <- function() {
      ini({
        ltlag <- log(0.2); lka <- log(1.15); lcl <- log(0.135); lv <- log(8)
        prop.err <- 0.15; add.err <- 0.6
        eta.tlag ~ 0.5; eta.ka ~ 0.5; eta.cl ~ 0.1; eta.v ~ 0.1
      })
      model({
        tlag <- exp(ltlag + eta.tlag); ka <- exp(lka + eta.ka)
        cl <- exp(lcl + eta.cl); v <- exp(lv + eta.v)
        d/dt(gut) <- -ka * gut
        d/dt(central) <- ka * gut - (cl / v) * central
        lag(gut) <- tlag
        cp <- central / v
        cp ~ prop(prop.err) + add(add.err)
      })
    }
    .d <- nlmixr2data::warfarin
    .d <- .d[.d$dvid == "cp", ]
    .rf <- function(sm) {
      .nlmixr(lag.mod, .d, est = "focei",
              foceiControl(sensMethod = sm, print = 0,
                           maxOuterIterations = 0L, covMethod = ""))
    }
    fwd <- .rf("forward")
    adj <- .rf("adjoint")
    ## adjoint runs (no crash) and its inner-optimized objective matches the
    ## analytic-jump forward objective
    expect_equal(adj$objf, fwd$objf, tolerance = 1e-2)
  })

  test_that("default sensMethod is 'default' across the nlm family", {
    expect_equal(nlmControl()$sensMethod, "default")
    expect_equal(foceiControl()$sensMethod, "default")
    expect_equal(nlminbControl()$sensMethod, "default")
    expect_equal(optimControl()$sensMethod, "default")
    expect_equal(n1qn1Control()$sensMethod, "default")
    expect_equal(lbfgsb3cControl()$sensMethod, "default")
  })

  test_that("'default' sensMethod defers to getOption('nlmixr2est.adjoint')", {
    .r <- function() {
      .ui <- .mkNlmUi(nlmControl())  # sensMethod = "default"
      nlmixr2est:::.nlmAdjointResolve(.ui)
    }
    withr::with_options(list(nlmixr2est.adjoint = "adjoint"), {
      expect_true(.r()$useAdjoint)          # theo 1-cmt: 2 states, adjoint applies
    })
    withr::with_options(list(nlmixr2est.adjoint = "forward"), {
      expect_false(.r()$useAdjoint)
    })
    ## the package default policy is "forward"
    withr::with_options(list(nlmixr2est.adjoint = NULL), {
      expect_false(.r()$useAdjoint)
    })
  })

  test_that("default nlm fit matches an explicit forward fit", {
    .d <- nlmixr2data::theo_sd
    ## the default policy is "forward", so the default resolves to forward
    dft <- .nlmixr(one.cmt, .d, est = "nlm", nlmControl(solveType = "grad", print = 0))
    fwd <- .nlmixr(one.cmt, .d, est = "nlm",
                   nlmControl(sensMethod = "forward", solveType = "grad", print = 0))
    expect_equal(dft$objf, fwd$objf, tolerance = 1e-4)
    expect_equal(setNames(dft$parFixedDf$Estimate, rownames(dft$parFixedDf)),
                 setNames(fwd$parFixedDf$Estimate, rownames(fwd$parFixedDf)),
                 tolerance = 1e-4)
  })

})
