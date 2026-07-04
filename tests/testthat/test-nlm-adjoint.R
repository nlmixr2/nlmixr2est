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

  test_that("sensMethod='auto' picks adjoint when thetas exceed states", {
    ## 4 estimated thetas, 2 ODE states -> adjoint
    expect_true(nlmixr2est:::.nlmAdjointResolve(
      .mkNlmUi(nlmControl(sensMethod = "auto")))$useAdjoint)
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
    adjAuto <- .rf("auto")

    expect_equal(adjExplicit$objf, fwd$objf, tolerance = 1e-4)
    expect_equal(adjStiff$objf, fwd$objf, tolerance = 1e-4)
    expect_equal(adjAuto$objf, fwd$objf, tolerance = 1e-4)

    expect_equal(.est(adjExplicit), .est(fwd), tolerance = 1e-4)
    expect_equal(.est(adjStiff), .est(fwd), tolerance = 1e-4)
    expect_equal(.est(adjAuto), .est(fwd), tolerance = 1e-4)
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

  test_that("default sensMethod is 'auto' across the nlm family", {
    expect_equal(nlmControl()$sensMethod, "auto")
    expect_equal(foceiControl()$sensMethod, "auto")
    expect_equal(nlminbControl()$sensMethod, "auto")
    expect_equal(optimControl()$sensMethod, "auto")
    expect_equal(n1qn1Control()$sensMethod, "auto")
    expect_equal(lbfgsb3cControl()$sensMethod, "auto")
  })

  test_that("default (auto) nlm fit matches an explicit forward fit", {
    .d <- nlmixr2data::theo_sd
    ## 4 thetas > 2 states, so the auto default resolves to adjoint
    dft <- .nlmixr(one.cmt, .d, est = "nlm", nlmControl(solveType = "grad", print = 0))
    fwd <- .nlmixr(one.cmt, .d, est = "nlm",
                   nlmControl(sensMethod = "forward", solveType = "grad", print = 0))
    expect_equal(dft$objf, fwd$objf, tolerance = 1e-4)
    expect_equal(setNames(dft$parFixedDf$Estimate, rownames(dft$parFixedDf)),
                 setNames(fwd$parFixedDf$Estimate, rownames(fwd$parFixedDf)),
                 tolerance = 1e-4)
  })

  test_that("focei adjoint inner sensitivities match the forward fit", {
    .d <- nlmixr2data::theo_sd
    ## 3 etas > 2 ODE states so the inner sens auto-resolves to adjoint
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
    .rf <- function(sm) {
      .nlmixr(mod, .d, est = "focei",
              foceiControl(sensMethod = sm, print = 0, maxOuterIterations = 15L))
    }
    fwd <- .rf("forward")
    adj <- .rf("adjoint")
    expect_equal(adj$objf, fwd$objf, tolerance = 1e-3)
    expect_equal(setNames(adj$parFixedDf$Estimate, rownames(adj$parFixedDf)),
                 setNames(fwd$parFixedDf$Estimate, rownames(fwd$parFixedDf)),
                 tolerance = 1e-3)
  })

})
