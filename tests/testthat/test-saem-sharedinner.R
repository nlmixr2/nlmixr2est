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
