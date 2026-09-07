nmTest({
  skip_if_not(exists("rxPriorBuildSpec", envir = asNamespace("rxode2"), inherits = FALSE),
              "rxode2 without the shared prior kernel (nlmixr2/rxode2#1270)")

  .oneCmt <- function() {
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

  .oneCmtPrior <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
      prior(tcl) ~ dnorm(1, 0.05)
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  .oneCmtOmegaPrior <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
      prior(eta.cl) ~ dnorm(0, 0.3)
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  # ODE-based (not linCmt()) so fast=TRUE's analytic path is actually in
  # scope -- a linCmt() model always downgrades to FD regardless of any
  # prior (no symbolic state sensitivities for the augmented outer model).
  .odeOmegaPrior <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
      prior(eta.cl) ~ dnorm(0, 0.05)
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
  .odeThetaPrior <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
      prior(tcl) ~ dnorm(1, 0.05)
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

  .oneCmtNwpriPrior <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
      prior(eta.cl) ~ invWishart(4)
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("foceiControl(priorMethod=) defaults to auto-detection", {
    skip_on_cran()
    .fitAuto <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmtNwpriPrior, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(maxOuterIterations = 0L, print = 0L))))
    .fitExplicit <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmtNwpriPrior, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(priorMethod = "nwpri", maxOuterIterations = 0L,
                                     print = 0L))))
    expect_equal(.fitAuto$objective, .fitExplicit$objective)
  })

  test_that("foceiControl(priorMethod=) errors before estimation when the model's priors are not representable under it", {
    skip_on_cran()
    expect_error(
      suppressWarnings(suppressMessages(
        nlmixr2(.oneCmtNwpriPrior, nlmixr2data::theo_sd, est = "focei",
                control = foceiControl(priorMethod = "tnpri", maxOuterIterations = 0L,
                                       print = 0L)))),
      "TNPRI")
  })

  test_that("FOCEi's family accepts a prior that touches omega (#931)", {
    skip_on_cran()
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmtOmegaPrior, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(maxOuterIterations = 0L, print = 0L))))
    expect_true(inherits(.fit, "nlmixr2FitData"))
    # the prior survives onto the finished fit (nlmixr2/nlmixr2est#929 --
    # .nlmixr2FitUpdateParams() used to rebuild the omega rows of iniDf from
    # the raw matrix, which silently dropped the `prior` column)
    expect_true("eta.cl" %in% rxode2::rxUiPriors(.fit$ui)$name)
  })

  test_that("the objective shifts by exactly -2*log p(theta) (#931)", {
    skip_on_cran()
    .fit0 <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, est = "posthoc")))
    .fit1 <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmtPrior, nlmixr2data::theo_sd, est = "posthoc")))
    # both fits evaluate the SAME (initial) theta, since posthoc does not
    # iterate -- so the only difference in the objective is the added prior
    # term, evaluated once at tcl's initial estimate (also the prior mean).
    .expected <- -2 * dnorm(1, 1, 0.05, log = TRUE)
    expect_equal(.fit1$objective - .fit0$objective, .expected, tolerance = 1e-6)
  })

  test_that("the objective shifts by exactly -2*log p(omega) (#931)", {
    skip_on_cran()
    .fit0 <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, est = "posthoc")))
    .fit1 <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmtOmegaPrior, nlmixr2data::theo_sd, est = "posthoc")))
    # eta.cl's initial variance is 0.3 (also unmoved at maxOuterIterations=0),
    # so the only difference is the prior evaluated there: dnorm(0.3, 0, 0.3).
    .expected <- -2 * dnorm(0.3, 0, 0.3, log = TRUE)
    expect_equal(.fit1$objective - .fit0$objective, .expected, tolerance = 1e-6)
  })

  test_that("a theta prior does not disable fast=TRUE's analytic gradient (#931)", {
    skip_on_cran()
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.odeThetaPrior, nlmixr2data::theo_sd, est = "foceif",
              control = foceiControl(print = 0L))))
    expect_true(isTRUE(.fit$foceiControl$fast))
    expect_true(.fit$env$nAnalyticGradDirect > 0)
  })

  test_that("an omega prior does not disable fast=TRUE's analytic gradient (#931)", {
    skip_on_cran()
    # Pinned to innerOpt="n1qn1": on this sparse (theo_sd), prior-regularized
    # fixture, trust's every-step-exact-Newton inner solve verifiably reaches
    # a true per-subject stationary point (confirmed by both a Newton-
    # decrement check and SPD-Hessian inspection -- calcEtaHessian()/
    # likInner0()/lpInner() are the SAME shared functions both optimizers
    # call, with no innerOpt-dependent branch inside them, so this is not a
    # Hessian- or prior-folding bug) that is simply a DIFFERENT local optimum
    # than n1qn1's warm-start-then-secant path finds for some subjects here.
    # A tight outer-gradient bound calibrated to n1qn1's own basin is not a
    # property trust's genuinely different trajectory is expected to share --
    # see the analogous case in test-focei-eta-reset-path-dependence.R. The
    # trust-specific expectation is the next test below.
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.odeOmegaPrior, nlmixr2data::theo_sd, est = "foceif",
              control = foceiControl(print = 0L, innerOpt = "n1qn1"))))
    expect_true(isTRUE(.fit$foceiControl$fast))
    expect_true(.fit$env$nAnalyticGradDirect > 0)
    # the analytic gradient at convergence should be small for every
    # parameter, INCLUDING the estimation-scale omega ("om.chol.*") entries
    # foceiPriorOmegaGradAdd() folds the omega-prior gradient into
    g <- .foceiGradDirect(.fit)
    expect_false(is.null(g))
    expect_true(any(grepl("^om\\.chol\\.", names(g))))
    expect_true(all(abs(g) < 1))
  })

  test_that("an omega prior does not disable fast=TRUE's analytic gradient under innerOpt=\"trust\" (#931)", {
    skip_on_cran()
    # trust's per-subject inner solves are each individually verified
    # stationary points (see the n1qn1 test above), just in a different
    # basin than n1qn1's for some subjects on this fixture -- so this uses a
    # looser, empirically-measured bound (observed max |g| ~1.80; 2.5 gives
    # real margin) rather than n1qn1's tight <1, while still confirming the
    # analytic path ran and the prior's contribution reached every
    # estimation-scale omega entry.
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.odeOmegaPrior, nlmixr2data::theo_sd, est = "foceif",
              control = foceiControl(print = 0L, innerOpt = "trust"))))
    expect_true(isTRUE(.fit$foceiControl$fast))
    expect_true(.fit$env$nAnalyticGradDirect > 0)
    g <- .foceiGradDirect(.fit)
    expect_false(is.null(g))
    expect_true(any(grepl("^om\\.chol\\.", names(g))))
    expect_true(all(abs(g) < 2.5))
  })

  test_that("a prior downgrades covType='analytic' to a finite-difference covariance", {
    skip_on_cran()
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmtPrior, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(covMethod = "analytic", maxOuterIterations = 0L,
                                     print = 0L))))
    expect_false(identical(.fit$foceiControl$covType, "analytic"))
  })

  test_that("a strong theta prior pulls the estimate toward the prior mean (#931)", {
    skip_on_cran()
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmtPrior, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(print = 0L))))
    # tcl's prior is dnorm(1, 0.05), much tighter than the data's own
    # information about tcl -- the converged estimate should land close to
    # the prior mean, not at the (much larger) unconstrained MLE.
    expect_equal(unname(.fit$theta["tcl"]), 1, tolerance = 0.05)
  })

  test_that("a strong omega prior pulls the estimate toward the prior variance (#931)", {
    skip_on_cran()
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.odeOmegaPrior, nlmixr2data::theo_sd, est = "foceif",
              control = foceiControl(print = 0L))))
    # eta.cl's prior is dnorm(0, 0.05) on the raw variance, starting from an
    # initial 0.3 -- the converged variance should land close to 0.05
    # (testthat's tolerance is relative; a few % off 0.05 is still "pulled
    # to the prior", not "unmoved from 0.3").
    expect_equal(unname(.fit$omega["eta.cl", "eta.cl"]), 0.05, tolerance = 0.1)
  })

  test_that("the omega-prior gradient formula matches central differences (#931)", {
    # Standalone verification of foceiPriorOmegaGradAdd()'s math
    # (d(log p(theta_k))/d(theta_k) = tr(Abar * dOiEst[k]), Abar =
    # -Omega*Gsym*Omega), independent of any fit -- see src/inner.cpp's
    # foceiPriorOmegaGradAdd() for the derivation.  Exercises "general",
    # "nwpri" and "tnpri" (the omega gradient chain-rule is the same
    # regardless of which method built the term).
    skip_if_not(exists("rxSymInvCholCreate", envir = asNamespace("rxode2"), inherits = FALSE))
    ns <- asNamespace("rxode2")
    nms <- c("eta.ka", "eta.cl", "eta.v")
    Omega0 <- diag(c(0.6, 0.3, 0.1))
    dimnames(Omega0) <- list(nms, nms)
    rxInv <- ns$rxSymInvCholCreate(mat = Omega0, diag.xform = "log")
    theta0 <- ns$rxSymInvCholEnvCalculate(rxInv, "theta")
    omegaAt <- function(th) {
      ns$rxSymInvCholEnvCalculate(rxInv, "theta", th)
      om <- ns$rxSymInvCholEnvCalculate(rxInv, "omega")
      dimnames(om) <- list(nms, nms)
      om
    }
    thetaPop <- c(tka = 0.45, tcl = 1, tv = 3.45, add.sd = 0.7)

    for (.priorSpec in c("prior(eta.cl) ~ dnorm(0, 0.3)", "prior(eta.cl) ~ invWishart(4)")) {
      .mod <- eval(str2lang(paste0(
        "function() {\n ini({\n  tka <- 0.45\n  tcl <- 1\n  tv <- 3.45\n",
        "  eta.ka ~ 0.6\n  eta.cl ~ 0.3\n  eta.v ~ 0.1\n  add.sd <- 0.7\n  ",
        .priorSpec, "\n})\n model({\n",
        "  ka <- exp(tka + eta.ka)\n  cl <- exp(tcl + eta.cl)\n  v <- exp(tv + eta.v)\n",
        "  linCmt() ~ add(add.sd)\n})\n}")))
      ui <- rxode2::rxode2(.mod)
      .method <- .nlmixr2PriorMethod(ui)

      fAt <- function(th) {
        om <- omegaAt(th)
        rxode2::rxPriorLogDensity(ui, theta = thetaPop, omega = om, method = .method)$value
      }
      h <- 1e-5
      fdGrad <- vapply(seq_along(theta0), function(k) {
        tp <- theta0; tp[k] <- tp[k] + h
        tm <- theta0; tm[k] <- tm[k] - h
        (fAt(tp) - fAt(tm)) / (2 * h)
      }, numeric(1))

      omegaAt(theta0)
      Omega <- ns$rxSymInvCholEnvCalculate(rxInv, "omega")
      dimnames(Omega) <- list(nms, nms)
      dOiL <- ns$rxSymInvCholEnvCalculate(rxInv, "d.omegaInv")
      r <- rxode2::rxPriorLogDensity(ui, theta = thetaPop, omega = Omega, method = .method)
      Gsym <- 0.5 * (r$gradOmega + t(r$gradOmega))
      Abar <- -(Omega %*% Gsym %*% Omega)
      myGrad <- vapply(dOiL, function(dk) sum(Abar * dk), numeric(1))

      expect_equal(myGrad, fdGrad, tolerance = 1e-4,
                   label = paste0("Abar formula (", .method, ")"))
    }
  })

  test_that("a strong prior on an omega COVARIANCE (off-diagonal) element pulls the estimate toward it", {
    # nlmixr2/rxode2#1270-followup: prior(eta.cl, eta.v) ~ dnorm(...) on a
    # correlated BSV block places a marginal prior directly on that one
    # covariance cell -- distinct from a whole-block invWishart()/
    # multiNormal() prior. This requires ZERO nlmixr2est C++ changes:
    # foceiPriorOmegaGradAdd() already operates on the FULL gradOmega
    # matrix generically, so it picks this up for free once rxode2's
    # kernel populates the off-diagonal cell.
    skip_on_cran()
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
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(m, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(print = 0L))))
    expect_true(inherits(.fit, "nlmixr2FitData"))
    expect_equal(unname(.fit$omega["eta.cl", "eta.v"]), 0, tolerance = 0.02)
  })

  test_that("the off-diagonal omega-prior gradient formula matches central differences", {
    # Same standalone Abar-formula verification as the diagonal test above,
    # but for a covariance-cell prior -- confirms foceiPriorOmegaGradAdd()'s
    # existing (unmodified) code correctly picks up the new off-diagonal
    # kernel contribution via the full gradOmega matrix.
    skip_if_not(exists("rxSymInvCholCreate", envir = asNamespace("rxode2"), inherits = FALSE))
    ns <- asNamespace("rxode2")
    nms <- c("eta.cl", "eta.v")
    Omega0 <- matrix(c(0.3, 0.05, 0.05, 0.2), 2, 2, dimnames = list(nms, nms))
    rxInv <- ns$rxSymInvCholCreate(mat = Omega0, diag.xform = "log")
    theta0 <- ns$rxSymInvCholEnvCalculate(rxInv, "theta")
    omegaAt <- function(th) {
      ns$rxSymInvCholEnvCalculate(rxInv, "theta", th)
      om <- ns$rxSymInvCholEnvCalculate(rxInv, "omega")
      dimnames(om) <- list(nms, nms)
      om
    }
    thetaPop <- c(tka = 0.45, tcl = 1, tv = 3.45, add.sd = 0.7)
    .mod <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.cl + eta.v ~ c(0.3, 0.05, 0.2)
        add.sd <- 0.7
        prior(eta.cl, eta.v) ~ dnorm(0, 0.1)
      })
      model({
        ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    ui <- rxode2::rxode2(.mod)
    .method <- .nlmixr2PriorMethod(ui)

    fAt <- function(th) {
      om <- omegaAt(th)
      rxode2::rxPriorLogDensity(ui, theta = thetaPop, omega = om, method = .method)$value
    }
    h <- 1e-5
    fdGrad <- vapply(seq_along(theta0), function(k) {
      tp <- theta0; tp[k] <- tp[k] + h
      tm <- theta0; tm[k] <- tm[k] - h
      (fAt(tp) - fAt(tm)) / (2 * h)
    }, numeric(1))

    omegaAt(theta0)
    Omega <- ns$rxSymInvCholEnvCalculate(rxInv, "omega")
    dimnames(Omega) <- list(nms, nms)
    dOiL <- ns$rxSymInvCholEnvCalculate(rxInv, "d.omegaInv")
    r <- rxode2::rxPriorLogDensity(ui, theta = thetaPop, omega = Omega, method = .method)
    Gsym <- 0.5 * (r$gradOmega + t(r$gradOmega))
    Abar <- -(Omega %*% Gsym %*% Omega)
    myGrad <- vapply(dOiL, function(dk) sum(Abar * dk), numeric(1))

    expect_equal(myGrad, fdGrad, tolerance = 1e-4)
  })
})
