nmTest({
  skip_if_not(exists("rxPriorBuildSpec", envir = asNamespace("rxode2"), inherits = FALSE),
              "rxode2 without the shared prior kernel (nlmixr2/rxode2#1270)")

  test_that("imp/impmap/qrpem declare 'general'-level prior support individually (#932)", {
    expect_identical(attr(nlmixr2Est.imp, "nlmixr2Priors"), "general")
    expect_identical(attr(nlmixr2Est.impmap, "nlmixr2Priors"), "general")
    expect_identical(attr(nlmixr2Est.qrpem, "nlmixr2Priors"), "general")
  })

  test_that("impmap accepts a prior that touches omega (#932)", {
    # the "theta"-level gate (a prior on omega is refused) is covered
    # generically for any method declaring that level in test-priors-assert.R;
    # impmap's own declared level is "general", so this just confirms the
    # fit runs and the prior survives onto the finished fit, mirroring
    # test-focei-prior.R's analogous check for FOCEi's family (#931).
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
    .fi <- suppressWarnings(suppressMessages(
      nlmixr2(m, nlmixr2data::theo_sd, est = "impmap",
              control = impmapControl(print = 0L, nIter = 1L))))
    expect_true(inherits(.fi, "nlmixr2FitCore"))
    expect_true("eta.cl" %in% rxode2::rxUiPriors(.fi$ui)$name)
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

  test_that("a strong invWishart prior pulls the Omega EM update toward its scale matrix (#932)", {
    # eta.cl's invWishart() term is CONJUGATE to the EM moment-average update
    # (impPriorOmegaCorrect()'s exact closed-form path, not the generic
    # Fisher-scoring one) -- rho=1000 dominates the data's own information
    # (12 subjects), so the converged variance should land close to the
    # term's own scale (the ini()-declared 0.3), not near the unconstrained
    # data-driven value.
    m <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3
        add.sd <- 0.7
        prior(eta.cl) ~ invWishart(1000)
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .fi <- suppressWarnings(suppressMessages(
      nlmixr2(m, nlmixr2data::theo_sd, est = "impmap",
              control = impmapControl(print = 0L, nIter = 30L, isample = 300L))))
    expect_true(inherits(.fi, "nlmixr2FitCore"))
    expect_equal(unname(.fi$omega["eta.cl", "eta.cl"]), 0.3, tolerance = 0.1)
  })

  test_that("a strong normal prior directly on omega pulls the EM update toward it (generic/tnpri path) (#932)", {
    # prior(eta.cl) ~ dnorm() on the raw omega value (not invWishart()) takes
    # impPriorOmegaCorrect()'s generic one-step Fisher-scoring branch, not the
    # invWishart conjugate closed form -- distinct code path from the test
    # above.  Starts from a much larger initial variance (0.6) than the
    # tight prior (0.05) it should be pulled toward.
    m <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.1; eta.cl ~ 0.6
        add.sd <- 0.7
        prior(eta.cl) ~ dnorm(0, 0.05)
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .fi <- suppressWarnings(suppressMessages(
      nlmixr2(m, nlmixr2data::theo_sd, est = "impmap",
              control = impmapControl(print = 0L, nIter = 30L, isample = 300L))))
    expect_true(inherits(.fi, "nlmixr2FitCore"))
    expect_equal(unname(.fi$omega["eta.cl", "eta.cl"]), 0.05, tolerance = 0.15)
    # Omega stays a valid (symmetric positive-definite) covariance under the
    # correction -- impSetOmega()'s PD floor/symmetrize guard degrades this
    # the same way an ordinary bad EM iterate would rather than producing NaN.
    expect_true(all(eigen(.fi$omega, symmetric = TRUE, only.values = TRUE)$values > 0))
  })

  test_that("a strong prior on an omega COVARIANCE (off-diagonal) element pulls the EM update toward it (generic path)", {
    # nlmixr2/rxode2#1270-followup: prior(eta.cl, eta.v) ~ dnorm(...) on a
    # correlated BSV block. impPriorOmegaCorrect()'s generic (non-invWishart)
    # branch operates on the full Omega/Gsym matrices already -- this
    # exercises that same code path, unmodified, for a covariance cell
    # rather than a diagonal one.
    skip_if_not(exists("rxPriorBuildSpec", envir = asNamespace("rxode2"), inherits = FALSE))
    m <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.cl + eta.v ~ c(0.3, 0.05, 0.2)
        add.sd <- 0.7
        prior(eta.cl, eta.v) ~ dnorm(0, 0.01)
      })
      model({
        ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .fi <- suppressWarnings(suppressMessages(
      nlmixr2(m, nlmixr2data::theo_sd, est = "impmap",
              control = impmapControl(print = 0L, nIter = 30L, isample = 300L))))
    expect_true(inherits(.fi, "nlmixr2FitCore"))
    expect_equal(unname(.fi$omega["eta.cl", "eta.v"]), 0, tolerance = 0.02)
    expect_true(all(eigen(.fi$omega, symmetric = TRUE, only.values = TRUE)$values > 0))
  })
})
