test_that("saemSharedResid_ (sharedInner='shared') computes f/r through likInner0", {
  skip_on_cran()

  # SAEM kernel unification: the shared FOCEi inner driver (likInner0) must
  # reproduce the per-observation prediction f and residual variance r that
  # SAEM otherwise computes with its own res_mod/arResk machinery.  Here we set
  # up the inner (as saemControl(sharedInner="shared") does) and check the
  # shared driver's f/r at eta=0 against an independent structural solve.
  m <- function() {
    ini({ tka <- log(1.5); tcl <- log(0.04); tv <- log(0.5); eta.cl ~ 0.1; add.sd <- 0.7 })
    model({
      ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  d <- nlmixr2data::theo_sd
  .N <- length(unique(d$ID))
  .zeroEta <- matrix(0.0, .N, 1L)
  .fc <- list(rxControl = rxode2::rxControl(), fastInnerIt = 100L, sumProd = FALSE,
              optExpression = TRUE, literalFix = FALSE, addProp = "combined2",
              eventSens = "jump", indTolRelax = TRUE, maxOdeRecalc = 5L,
              odeRecalcFactor = 10^0.5)
  .env <- .fsaemInnerSetup(m(), d, .zeroEta, .fc)
  .res <- saemSharedResid_(.zeroEta)
  .fsaemInnerFree()

  # independent structural solve at eta=0 (plain rxode2, no BSV)
  .ms <- rxode2::rxode2({
    ka <- exp(log(1.5)); cl <- exp(log(0.04)); v <- exp(log(0.5))
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
  })
  .fIndep <- as.data.frame(rxode2::rxSolve(.ms, d, addDosing = FALSE))$cp

  expect_equal(length(.res$f), length(.fIndep))
  # f from the shared inner matches the independent solve (ODE-tolerance level)
  expect_lt(max(abs(.res$f - .fIndep)), 0.05)
  # additive error: r is constant at add.sd^2 across all observations
  expect_true(all(abs(.res$r - 0.7^2) < 1e-4))
  expect_lt(diff(range(.res$r)), 1e-6)
  expect_true(is.finite(.res$objf))
})

test_that("shared inner driver reproduces a converged SAEM fit's IPRED / r", {
  skip_on_cran()

  # Full-fit equivalence: at the CONVERGED theta/omega/etas, the shared driver's
  # per-observation f must reproduce the fit's individual predictions (IPRED),
  # and r the additive variance -- i.e. saemSharedResid_ is a drop-in for SAEM's
  # residual/prediction computation.
  m <- function() {
    ini({ tka <- log(1.5); tcl <- log(0.04); tv <- log(0.5); eta.cl ~ 0.1; add.sd <- 0.7 })
    model({
      ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  d <- nlmixr2data::theo_sd
  .f <- suppressWarnings(nlmixr2(m, d, est = "saem",
    control = saemControl(nBurn = 100, nEm = 50, print = 0)))

  .th <- .f$ui$iniDf
  .th <- .th[!is.na(.th$ntheta), ]
  .th <- .th[order(.th$ntheta), ]
  .thetaFull <- as.numeric(.f$theta[.th$name])
  .omega <- diag(.f$omega)
  .condEta <- as.matrix(.f$eta[order(.f$eta$ID), "eta.cl", drop = FALSE])

  .fc <- list(rxControl = rxode2::rxControl(), fastInnerIt = 100L, sumProd = FALSE,
              optExpression = TRUE, literalFix = FALSE, addProp = "combined2",
              eventSens = "jump", indTolRelax = TRUE, maxOdeRecalc = 5L,
              odeRecalcFactor = 10^0.5)
  .env <- .fsaemInnerSetup(m(), d, .condEta, .fc)
  .fsaemInnerUpdate(.env, .thetaFull, .omega, .condEta)
  .res <- saemSharedResid_(.condEta)
  .fsaemInnerFree()

  expect_equal(length(.res$f), length(.f$IPRED))
  # the shared driver's f reproduces the fit's IPRED (ODE-tolerance level)
  expect_lt(max(abs(.res$f - .f$IPRED)), 0.05)
  # r is the additive residual variance (add.sd^2) at the converged add.sd
  .addSd <- as.numeric(.f$theta[.th$name])[.th$name == "add.sd"]
  expect_true(all(abs(.res$r - .addSd^2) < 1e-3))
})

test_that("shared inner driver handles a PROPORTIONAL error model (r = (prop.sd*f)^2)", {
  skip_on_cran()

  # The driver reads r from rxInner's rx_r_, so it is error-model-agnostic;
  # confirm proportional error reproduces r = (prop.sd*f)^2.
  m <- function() {
    ini({ tka <- log(1.5); tcl <- log(0.04); tv <- log(0.5); eta.cl ~ 0.1; prop.sd <- 0.2 })
    model({
      ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ prop(prop.sd)
    })
  }
  d <- nlmixr2data::theo_sd
  .f <- suppressWarnings(nlmixr2(m, d, est = "saem",
    control = saemControl(nBurn = 100, nEm = 50, print = 0)))
  .th <- .f$ui$iniDf; .th <- .th[!is.na(.th$ntheta), ]; .th <- .th[order(.th$ntheta), ]
  .thetaFull <- as.numeric(.f$theta[.th$name]); .omega <- diag(.f$omega)
  .condEta <- as.matrix(.f$eta[order(.f$eta$ID), "eta.cl", drop = FALSE])
  .fc <- list(rxControl = rxode2::rxControl(), fastInnerIt = 100L, sumProd = FALSE,
              optExpression = TRUE, literalFix = FALSE, addProp = "combined2",
              eventSens = "jump", indTolRelax = TRUE, maxOdeRecalc = 5L,
              odeRecalcFactor = 10^0.5)
  .env <- .fsaemInnerSetup(m(), d, .condEta, .fc)
  .fsaemInnerUpdate(.env, .thetaFull, .omega, .condEta)
  .res <- saemSharedResid_(.condEta)
  .fsaemInnerFree()

  expect_lt(max(abs(.res$f - .f$IPRED)), 0.05)
  .propSd <- as.numeric(.f$theta[.th$name])[.th$name == "prop.sd"]
  .rExpect <- (.propSd * .res$f)^2
  .keep <- .res$f > 1   # ignore near-zero predictions
  expect_lt(max(abs(.res$r[.keep] - .rExpect[.keep]) / .rExpect[.keep]), 1e-3)
})
