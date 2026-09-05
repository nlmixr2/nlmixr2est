nmTest({
  test_that("foceiControl(innerOpt=) trust mapping", {
    expect_equal(foceiControl()$innerOpt, 4L)
    expect_equal(foceiControl(innerOpt = "n1qn1")$innerOpt, 1L)
    expect_equal(foceiControl(innerOpt = "BFGS")$innerOpt, 2L)
    expect_equal(foceiControl(innerOpt = "trust")$innerOpt, 3L)
    expect_equal(foceiControl(innerOpt = "auto")$innerOpt, 4L)
    expect_equal(foceiControl(innerOpt = 3L)$innerOpt, 3L)
    expect_error(foceiControl(innerOpt = "bogus"))
    expect_error(foceiControl(innerOpt = 5L))

    expect_equal(foceiControl()$trustConf, 0.975)
    expect_equal(foceiControl(trustConf = 0.9)$trustConf, 0.9)
    expect_null(foceiControl()$trustRinit)
    expect_null(foceiControl()$trustRmax)
    expect_equal(foceiControl(trustRmax = 2)$trustRmax, 2)
    expect_error(foceiControl(trustRinit = 5, trustRmax = 1))

    # trustConf must be strictly inside (0, 1): qchisq(0, df)==0 (zero
    # trust-region radius) and qchisq(1, df)==Inf (unbounded radius) are both
    # degenerate.
    expect_error(foceiControl(trustConf = 1))
    expect_error(foceiControl(trustConf = 0))

    # trustRinit/trustRmax==0 is the same degenerate case (a zero-radius trust
    # region that can never step).
    expect_error(foceiControl(trustRinit = 0))
    expect_error(foceiControl(trustRmax = 0))

    # trustFterm/trustMterm default to 10^(-sigdig-2), two orders tighter than
    # the plain 10^(-sigdig): the inner solve is the function the outer problem
    # differentiates, so this tolerance is the objective's noise floor and a
    # finite-difference outer gradient cannot resolve a step below it.  They
    # stay independent of epsilon (n1qn1's own, unrelated tolerance).
    expect_equal(foceiControl(sigdig = 3)$trustFterm, 1e-5)
    expect_equal(foceiControl(sigdig = 3)$trustMterm, 1e-5)
    expect_equal(foceiControl(sigdig = 4)$trustFterm, 1e-6)
    expect_equal(foceiControl(sigdig = 3, epsilon = 1e-2)$trustFterm, 1e-5)
    expect_equal(foceiControl(sigdig = 3, epsilon = 1e-2)$trustMterm, 1e-5)
    expect_equal(foceiControl(trustFterm = 0.5)$trustFterm, 0.5)
    expect_equal(foceiControl(trustMterm = 0.25)$trustMterm, 0.25)
    expect_error(foceiControl(trustFterm = 0))
    expect_error(foceiControl(trustMterm = 0))
    expect_error(foceiControl(trustFterm = -1))
    expect_error(foceiControl(trustMterm = -1))

    .ctl <- foceiControl(innerOpt = "trust", trustConf = 0.9,
                         trustFterm = 0.5, trustMterm = 0.25)
    expect_equal(do.call(foceiControl, .ctl)$innerOpt, 3L)
    expect_equal(do.call(foceiControl, .ctl)$trustConf, 0.9)
    expect_equal(do.call(foceiControl, .ctl)$trustFterm, 0.5)
    expect_equal(do.call(foceiControl, .ctl)$trustMterm, 0.25)
  })

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

  .fitTrustCmp <- function(innerOpt) {
    suppressWarnings(suppressMessages(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(innerOpt = innerOpt, maxOuterIterations = 20,
                                      covMethod = "", calcTables = FALSE, print = 0))
    ))
  }

  test_that("innerOpt='trust' converges close to n1qn1", {
    skip_on_cran()
    .f1 <- .fitTrustCmp("n1qn1")
    .n1 <- .nTrustInner()
    .f2 <- .fitTrustCmp("trust")
    .n2 <- .nTrustInner()

    expect_true(is.finite(.f1$objf))
    expect_true(is.finite(.f2$objf))
    expect_equal(.f2$objf, .f1$objf, tolerance = 1e-1)
    expect_equal(as.data.frame(.f1$eta), as.data.frame(.f2$eta), tolerance = 5e-2)

    # Positive evidence the trust path actually ran -- not a silent fallback to
    # n1qn1, the failure mode #927's innerOpt="BFGS" had (numeric agreement alone
    # would not catch that).
    expect_equal(.n1, 0L)
    expect_true(.n2 > 0L)
  })

  test_that("innerOpt='trust' clamps trustRinit to a smaller derived trustRmax", {
    skip_on_cran()
    # trustRmax's default depends on neta (only known in C++), so R can only
    # reject trustRinit > trustRmax when BOTH are given explicitly. Leaving
    # trustRmax NULL with a large explicit trustRinit reaches that same
    # geometry violation (trustRinit > trustRmax) via the DERIVED default
    # instead -- the C++ side clamps trustRinit down rather than starting the
    # trust region already past its own cap. This just has to not error/hang.
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(innerOpt = "trust", trustRinit = 100,
                                      maxOuterIterations = 20,
                                      covMethod = "", calcTables = FALSE, print = 0))
    ))
    expect_true(is.finite(.fit$objf))
  })

  test_that("innerOpt='BFGS' actually falls back to n1qn1 (not just the R-level mapping)", {
    skip_on_cran()
    # #927: innerOpt="BFGS" is accepted but unimplemented in C++ (lbfgsb3C is not
    # reentrant under this OpenMP loop, see src/inner.cpp). Run a real fit, not just
    # check foceiControl()$innerOpt, so a future C++ change that actually wires
    # innerOpt==2 into the trust/lbfgsb3C path gets caught here too.
    .f1 <- .fitTrustCmp("n1qn1")
    .fB <- .fitTrustCmp("BFGS")
    .nB <- .nTrustInner()

    expect_equal(.fB$objf, .f1$objf, tolerance = 1e-8)
    expect_equal(as.data.frame(.fB$eta), as.data.frame(.f1$eta), tolerance = 1e-8)
    expect_equal(.nB, 0L)
  })

  test_that("innerOpt='auto' picks n1qn1 for a generalized likelihood and trust otherwise", {
    skip_on_cran()
    # "auto" resolves in C++ (foceiSetup_), the first point where needOptimHess
    # is known, so assert on .nTrustInner() -- the R control still reports 4L.
    .fN <- .fitTrustCmp("auto")
    expect_true(.nTrustInner() > 0L)

    .llCmt <- .oneCmt |> model(linCmt() ~ add(add.sd) + dnorm())
    .fLL <- suppressWarnings(suppressMessages(
      nlmixr2(.llCmt, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(innerOpt = "auto", maxOuterIterations = 20,
                                      covMethod = "", calcTables = FALSE, print = 0))
    ))
    expect_equal(.nTrustInner(), 0L)
    expect_true(is.finite(.fLL$objf))
  })

  test_that("innerOpt='trust' restart cascade completes when the inner solve doesn't converge", {
    skip_on_cran()
    # A tiny maxInnerIterations forces trust_solve_c() to hit iterlim before
    # converging (converged=FALSE, no NA involved), exercising the !converged
    # nudge-restart cascade in the trustInner branch of innerOpt1()
    # (src/inner.cpp) end to end: this only proves the cascade runs and still
    # returns a usable fit, NOT specifically that fInd->badSolve is reset per
    # attempt (an outside review caught that trustSolveAt() must do this --
    # without it, an NA during the FIRST attempt permanently poisons every
    # later nudge via trustInnerObjfun's own badSolve guard). Forcing a real
    # NA reliably and portably (vs. this iterlim path) would need a
    # deliberately pathological model; the badSolve reset itself is a
    # one-line, easily re-verified-by-reading fix instead.
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(innerOpt = "trust", maxOuterIterations = 5,
                                      maxInnerIterations = 2,
                                      covMethod = "", calcTables = FALSE, print = 0))
    ))
    expect_true(is.finite(.fit$objf))
    expect_true(.nTrustInner() > 0L)
  })

  # A non-normal general-ll() endpoint: calcEtaHessian()'s needOptimHess
  # branch is a genuinely different code path from the normal-endpoint
  # analytic Gauss-Newton branch the rest of this file's models exercise
  # (see calcEtaHessian(), src/inner.cpp) -- innerOpt="trust" must not be
  # implicitly Gaussian-family-only.
  .poisMod <- function() {
    ini({
      te0 <- log(5)
      eta.e0 ~ 0.3
    })
    model({
      e0 <- exp(te0 + eta.e0)
      effect <- e0
      effect ~ dpois(effect)
    })
  }
  .poisData <- data.frame(ID = rep(1:20, each = 3), TIME = rep(c(0, 1, 2), 20))
  .poisData$DV <- c(5L, 4L, 6L, 3L, 5L, 4L, 6L, 5L, 7L, 4L, 5L, 5L,
                     6L, 6L, 5L, 4L, 3L, 5L, 5L, 6L, 4L, 7L, 5L, 6L,
                     4L, 5L, 5L, 3L, 4L, 6L, 6L, 5L, 4L, 5L, 6L, 5L,
                     4L, 4L, 6L, 5L, 7L, 5L, 6L, 5L, 4L, 3L, 5L, 6L,
                     5L, 5L, 6L, 4L, 6L, 5L, 5L, 4L, 5L, 6L, 4L, 5L)

  test_that("innerOpt='trust' converges close to n1qn1 on a non-normal (dpois) endpoint", {
    skip_on_cran()
    .f1 <- .nlmixr(.poisMod, .poisData, est = "focei",
                   control = foceiControl(innerOpt = "n1qn1", print = 0L))
    .n1 <- .nTrustInner()
    .f2 <- .nlmixr(.poisMod, .poisData, est = "focei",
                   control = foceiControl(innerOpt = "trust", print = 0L))
    .n2 <- .nTrustInner()

    expect_true(is.finite(.f1$objf))
    expect_true(is.finite(.f2$objf))
    expect_equal(.f2$objf, .f1$objf, tolerance = 1e-1)
    expect_equal(unname(.f1$theta), unname(.f2$theta), tolerance = 5e-2)

    expect_equal(.n1, 0L)
    expect_true(.n2 > 0L)
  })

  test_that("a fit reports whether its trust inner solves converged (#1044)", {
    skip_on_cran()
    # op_focei.nTrustInner counted trust_solve_c() CALLS, so a fit whose inner
    # solves all failed was indistinguishable, from the fit object, from one
    # where they all converged -- which is what made #1040's set B impossible to
    # diagnose.  The outcome breakdown is what answers that.
    .ok <- .fitTrustCmp("trust")
    .cnt <- .ok$env$nTrustInner
    expect_true(is.integer(.cnt))
    expect_equal(sort(names(.cnt)),
                 sort(c("calls", "error", "notConverged", "solverFail",
                        "newtonGate", "warmRetry", "radiusRetry", "nudge",
                        "failed")))
    expect_gt(.cnt[["calls"]], 0L)
    # This fit converges cleanly, so nothing below "calls" fires.  That is the
    # half of the diagnostic that has to stay quiet or it says nothing.
    expect_equal(.cnt[["failed"]], 0L)
    expect_equal(.cnt[["notConverged"]], 0L)
    expect_equal(.cnt[["nudge"]], 0L)

    # n1qn1 has no trust solves at all, so the entry is absent rather than zero.
    expect_null(.fitTrustCmp("n1qn1")$env$nTrustInner)

    # maxInnerIterations=2 makes trust_solve_c() hit iterlim, so it reports the
    # non-convergence itself (solverFail) rather than the Newton-decrement gate
    # withdrawing it, and some subjects exhaust the whole cascade.
    .bad <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(innerOpt = "trust", maxOuterIterations = 5,
                                     maxInnerIterations = 2, covMethod = "",
                                     calcTables = FALSE, print = 0))))
    .bc <- .bad$env$nTrustInner
    expect_gt(.bc[["solverFail"]], 0L)
    # An attempt ends non-converged in exactly one of three ways, so the
    # breakdown has to add up or it is not a breakdown.
    expect_equal(.bc[["notConverged"]],
                 .bc[["error"]] + .bc[["solverFail"]] + .bc[["newtonGate"]])
    expect_gt(.bc[["failed"]], 0L)
    expect_true(is.finite(.bad$objf))
  })

  test_that("a failed inner attempt cannot win the marginal re-rank (#1044)", {
    skip_on_cran()
    # The restart candidates the marginal re-rank chooses from carried no record
    # of whether the attempt that produced them succeeded, so a failed attempt's
    # eta could be selected over a converged one -- its Laplace log|H| is
    # measured at a point the inner objective never descended to.  "dropped" is
    # the count that shows the rule fires; equivalence of objectives cannot
    # distinguish "dropped a bad candidate" from "never had one".
    .bad <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(innerOpt = "trust", maxOuterIterations = 5,
                                     maxInnerIterations = 2, covMethod = "",
                                     calcTables = FALSE, print = 0))))
    .rr <- .bad$env$nInnerRerank
    expect_equal(sort(names(.rr)),
                 sort(c("ranked", "flipped", "noGood", "dropped")))
    expect_gt(.rr[["dropped"]], 0L)
    # Solves where EVERY candidate failed still have to report something, so the
    # rule falls back to them rather than returning nothing.
    expect_gt(.rr[["noGood"]], 0L)
    expect_true(is.finite(.bad$objf))

    # A fit whose inner solves all converge drops nothing, so the rule costs the
    # ordinary path neither a candidate nor an extra evaluation.
    .ok <- .fitTrustCmp("trust")
    expect_equal(.ok$env$nInnerRerank[["dropped"]], 0L)
    expect_equal(.ok$env$nInnerRerank[["noGood"]], 0L)
  })

  # A model with curvature the inner problem can get lost in -- two etas with a
  # FIXED unit omega entering through an inverse CDF, alongside an estimated
  # block (the #1040 model).  Displacing every theta by `d` walks the inner
  # solves from "all converge" into "most trip the Newton-decrement gate", which
  # is the regime #1044 is about and the only one the retry counts fire in.
  .gateBase <- function() {
    ini({
      tcl <- log(4)
      tv1 <- log(30)
      tq <- log(4)
      tv2 <- log(40)
      rxz.eta.cl ~ fix(1)
      rxz.eta.v1 ~ fix(1)
      eta.q + eta.v2 ~ c(0.0305,
                         0.0107, 0.0285)
      prop.sd <- 0.1
    })
    model({
      cl <- exp(tcl + 0.3 * logit(pnorm(rxz.eta.cl)))
      v1 <- exp(tv1 + 0.3 * logit(pnorm(rxz.eta.v1)))
      q <- exp(tq + eta.q)
      v2 <- exp(tv2 + eta.v2)
      linCmt() ~ prop(prop.sd)
    })
  }

  .gateMod <- function(d) {
    suppressWarnings(suppressMessages(
      rxode2::ini(.gateBase, tcl = log(4) + d, tv1 = log(30) - d,
                  tq = log(4) + 1.5 * d, tv2 = log(40) + d)))
  }

  .gateData <- function() {
    withr::with_seed(42, {
      .ev <- rxode2::et(amt = 200, ii = 24, addl = 2)
      .ev <- rxode2::et(.ev, seq(0.5, 96, length.out = 12))
      .ev <- rxode2::et(.ev, id = 1:60)
      .d <- suppressWarnings(suppressMessages(
        as.data.frame(rxode2::rxSolve(.gateBase, .ev, addDosing = TRUE))))
    })
    .dat <- data.frame(ID = .d$id, TIME = .d$time, DV = .d$sim, EVID = .d$evid,
                       AMT = ifelse(is.na(.d$amt), 0, .d$amt))
    .dat$DV[.dat$EVID != 0] <- NA
    .dat
  }

  test_that("the trust retry cascade's stages are each reachable (#1044)", {
    skip_on_cran()
    .dat <- .gateData()
    .gateFit <- function(d) {
      suppressWarnings(suppressMessages(
        nlmixr2(.gateMod(d), .dat, "focei",
                foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L,
                             maxInnerIterations = 5000L, calcTables = FALSE,
                             innerOpt = "trust"))))
    }

    # Each retry stage needs its own count > 0 or it could be dead code and
    # every other assertion here would still pass.  Only ">0" is asserted: the
    # exact numbers move with rxode2's solver.
    .g3 <- .gateFit(3)
    .c3 <- .g3$env$nTrustInner
    # trust_solve_c() reports convergence and the Newton decrement says
    # otherwise -- this is what drives the retries, not solverFail.
    expect_gt(.c3[["newtonGate"]], 0L)
    # The in-place re-solve, taken when the Newton step still fits the radius.
    expect_gt(.c3[["warmRetry"]], 0L)
    expect_gt(.c3[["nudge"]], 0L)
    # ... and subjects where none of it worked, which is the count that says a
    # fit is in trouble.
    expect_gt(.c3[["failed"]], 0L)
    expect_true(is.finite(.g3$objf))

    # Displaced further, the Newton step outgrows the radius often enough to
    # exercise the escalation branch instead.
    .c4 <- .gateFit(4)$env$nTrustInner
    expect_gt(.c4[["radiusRetry"]], 0L)
  })

  test_that("mceta's mixed candidate set drops the failed attempts (#1044)", {
    skip_on_cran()
    # The re-rank only has something to choose between when the restarts
    # produce several candidates, and mceta >= 1 is what makes that ordinary.
    # maxInnerIterations=2 makes some of those attempts genuinely fail, so the
    # candidate set is mixed -- which is the case the rule exists for.
    .mc <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(innerOpt = "trust", maxOuterIterations = 5,
                                     maxInnerIterations = 2, mceta = 5L,
                                     covMethod = "", calcTables = FALSE,
                                     print = 0))))
    expect_gt(.mc$env$nTrustInner[["solverFail"]], 0L)
    expect_gt(.mc$env$nInnerRerank[["dropped"]], 0L)
    expect_gt(.mc$env$nInnerRerank[["ranked"]], 0L)
    expect_true(is.finite(.mc$objf))

    # The same fit with converging inner solves drops nothing, so mceta does not
    # pay for the rule when there is nothing to drop.
    .ok <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(innerOpt = "trust", maxOuterIterations = 5,
                                     mceta = 5L, covMethod = "",
                                     calcTables = FALSE, print = 0))))
    expect_gt(.ok$env$nInnerRerank[["ranked"]], 0L)
    expect_equal(.ok$env$nInnerRerank[["dropped"]], 0L)
  })

  test_that("innerOpt='trust' is thread-safe (cores>=2 matches serial)", {
    skip_on_cran()
    .old <- rxode2::getRxThreads()
    on.exit(rxode2::setRxThreads(.old), add = TRUE)

    rxode2::setRxThreads(1L)
    .f1 <- .fitTrustCmp("trust")

    rxode2::setRxThreads(2L)
    .f2 <- .fitTrustCmp("trust")

    expect_equal(.f2$objf, .f1$objf, tolerance = 1e-8)
    expect_equal(as.data.frame(.f1$eta), as.data.frame(.f2$eta), tolerance = 1e-8)
  })
})
