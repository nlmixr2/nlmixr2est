nmTest({
  skip_if_not(exists("rxPriorBuildSpec", envir = asNamespace("rxode2"), inherits = FALSE),
              "rxode2 without the shared prior kernel (nlmixr2/rxode2#1270)")

  test_that("imp/impmap/qrpem declare 'theta'-level prior support individually (#932)", {
    expect_identical(attr(nlmixr2Est.imp, "nlmixr2Priors"), "theta")
    expect_identical(attr(nlmixr2Est.impmap, "nlmixr2Priors"), "theta")
    expect_identical(attr(nlmixr2Est.qrpem, "nlmixr2Priors"), "theta")
  })

  test_that("impmap refuses a prior touching omega at the 'theta' level (#932)", {
    m <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3
        add.sd <- 0.7
        prior(eta.cl) ~ dnorm(0, 0.3)
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    expect_error(
      suppressWarnings(suppressMessages(
        nlmixr2(m, nlmixr2data::theo_sd, est = "impmap",
                control = impmapControl(print = 0L, nIter = 1L)))),
      "eta.cl")
  })

  test_that("a strong prior on a non-mu structural theta pulls the M-step Newton estimate toward it (#932)", {
    # add.sd is a non-mu structural (residual-error) theta, estimated exclusively
    # by impUpdateStructThetas()'s Newton step on the theta-sensitivity score/
    # Hessian (impmapControl "M6") -- the finalize MAP-pass recompute never
    # touches it, so this exercises impPriorStructThetaCorrect() specifically.
    m <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.cl ~ 0.1
        add.sd <- 0.7
        prior(add.sd) ~ dnorm(2, 0.02)
      })
      model({
        ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d / dt(depot) <- -ka * depot
        d / dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp ~ add(add.sd)
      })
    }
    .fi <- suppressWarnings(suppressMessages(
      nlmixr2(m, nlmixr2data::theo_sd, est = "impmap",
              control = impmapControl(print = 0L, nIter = 30L, isample = 300L))))
    expect_true(inherits(.fi, "nlmixr2FitCore"))
    expect_equal(unname(fixef(.fi)["add.sd"]), 2, tolerance = 0.15)
  })

  test_that("a strong prior on a plain mu-intercept theta pulls the mean-shift EM estimate toward it (#932)", {
    # tcl is a plain (no-covariate) mu-referenced intercept, handled exclusively
    # by impMuInterceptStep()'s mean-shift EM update -- distinct from both the
    # struct-theta Newton step above and the covariate regression below.
    m <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- fix(3.45)
        eta.ka ~ 0.6; eta.cl ~ 0.3
        add.sd <- fix(0.7)
        prior(tcl) ~ dnorm(3, 0.02)
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .fi <- suppressWarnings(suppressMessages(
      nlmixr2(m, nlmixr2data::theo_sd, est = "imp",
              control = impControl(print = 0L, nIter = 40L, isample = 300L))))
    expect_true(inherits(.fi, "nlmixr2FitCore"))
    expect_equal(unname(fixef(.fi)["tcl"]), 3, tolerance = 0.1)
  })

  test_that("a strong prior on a mu-referenced covariate coefficient pulls updateMuGroups() toward it (#932)", {
    # cl.wt is a mu-referenced covariate coefficient, handled exclusively by
    # updateMuGroups()'s (weighted) normal-equations regression -- a real WT
    # effect is in the data, so without a prior this lands well away from 0
    # (see test-impmap.R's "M4: mu-referenced covariate" for the unconstrained
    # FOCEI estimate); a tight prior(cl.wt) ~ dnorm(0, ...) should pull it back.
    m <- function() {
      ini({
        tka <- 0.45; tcl <- 1; cl.wt <- 0.75; tv <- fix(3.45)
        eta.ka ~ 0.6; eta.cl ~ 0.3
        add.sd <- fix(0.7)
        prior(cl.wt) ~ dnorm(0, 0.001)
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + cl.wt * log(WT / 70) + eta.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .fi <- suppressWarnings(suppressMessages(
      nlmixr2(m, nlmixr2data::theo_sd, est = "impmap",
              control = impmapControl(print = 0L, nIter = 40L, isample = 300L))))
    expect_true(inherits(.fi, "nlmixr2FitCore"))
    expect_true(.fi$env$impMuGroupN >= 1L)
    expect_equal(unname(fixef(.fi)["cl.wt"]), 0, tolerance = 0.05)
  })
})
