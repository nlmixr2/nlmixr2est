nmTest({
  test_that("trustControl() defaults and validation", {
    expect_equal(trustControl()$fterm, 1e-5)
    expect_equal(trustControl()$mterm, 1e-5)
    expect_equal(trustControl(sigdig = 4)$fterm, 1e-6)
    expect_equal(trustControl(fterm = 0.5)$mterm, 0.5)
    expect_equal(trustControl(fterm = 0.5, mterm = 0.1)$mterm, 0.1)
    expect_null(trustControl()$rinit)
    expect_null(trustControl()$rmax)
    expect_equal(trustControl(rmax = 2)$rmax, 2)

    expect_error(trustControl(rinit = 0))
    expect_error(trustControl(rmax = 0))
    expect_error(trustControl(rinit = 5, rmax = 1))

    .ctl <- trustControl(rinit = 0.1, rmax = 1)
    expect_equal(do.call(trustControl, .ctl)$rinit, 0.1)
    expect_equal(do.call(trustControl, .ctl)$rmax, 1)

    expect_equal(trustControl()$hessianMethod, 3L)
    expect_equal(trustControl(hessianMethod = "fd")$hessianMethod, 1L)
    expect_equal(trustControl(hessianMethod = "bfgs")$hessianMethod, 2L)
    expect_equal(trustControl(hessianMethod = "sr1")$hessianMethod, 3L)
    expect_equal(trustControl(hessianMethod = "bofill")$hessianMethod, 4L)
    expect_equal(trustControl(hessianMethod = 3L)$hessianMethod, 3L)
    expect_error(trustControl(hessianMethod = "not-a-method"))
    expect_identical(trustControl()$covMethod, "r")
    expect_identical(trustControl(covMethod = "trust")$covMethod, "trust")
  })

  test_that("trust registers in the nlm-family method listing", {
    .types <- nlmixr2AllEstType()
    expect_true("trust" %in% .types$est)
    expect_equal(.types$type[.types$est == "trust"], "Optimizer (NLM family)")
  })

  .oneCmt <- function() {
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
      d / dt(depot) <- -ka * depot
      d / dt(centr) <- ka * depot - cl / v * centr
      cp <- centr / v
      cp ~ add(add.sd)
    })
  }

  test_that("est='trust' converges close to bobyqa", {
    skip_on_cran()
    .fB <- .nlmixr(
      .oneCmt,
      nlmixr2data::theo_sd,
      est = "bobyqa",
      control = bobyqaControl(print = 0L, calcTables = FALSE)
    )
    .fT <- .nlmixr(.oneCmt, nlmixr2data::theo_sd, est = "trust", control = trustControl(print = 0L, calcTables = FALSE))
    .nT <- .nTrustOuter()

    expect_true(is.finite(.fB$objective))
    expect_true(is.finite(.fT$objective))
    expect_equal(.fT$objective, .fB$objective, tolerance = 1e-2)
    expect_equal(unname(.fT$theta), unname(.fB$theta), tolerance = 1e-2)

    # Positive evidence the C++-resident loop actually ran -- not just that
    # the numbers happen to agree (mirrors test-focei-trust-inner.R's own
    # .nTrustInner() convention for the same reason).
    expect_true(.nT > 0L)
  })

  test_that("est='trust' converges close to bobyqa for every hessianMethod", {
    skip_on_cran()
    .fB <- .nlmixr(
      .oneCmt,
      nlmixr2data::theo_sd,
      est = "bobyqa",
      control = bobyqaControl(print = 0L, calcTables = FALSE)
    )
    for (.hm in c("fd", "bfgs", "sr1", "bofill")) {
      .fT <- .nlmixr(
        .oneCmt,
        nlmixr2data::theo_sd,
        est = "trust",
        control = trustControl(print = 0L, calcTables = FALSE, hessianMethod = .hm)
      )
      expect_true(is.finite(.fT$objective), info = .hm)
      expect_equal(.fT$objective, .fB$objective, tolerance = 1e-2, info = .hm)
      expect_equal(unname(.fT$theta), unname(.fB$theta), tolerance = 1e-2, info = .hm)
      expect_true(.nTrustOuter() > 0L, info = .hm)
    }
  })

  .oneCmtBounded <- function() {
    ini({
      tka <- 0.45
      tcl <- log(c(0, 2.7, 100))
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("est='trust' converges close to bobyqa on a bounded/upper_exp-transformed theta (#994)", {
    skip_on_cran()
    # tcl's internal coordinate has a near-zero starting gradient at this
    # model's default initial estimate, which used to make
    # nlmGetScaleC()'s unguarded derivative-based scaleC formula blow up to
    # scaleCmax (src/nlm.cpp) and permanently corrupt every later scaled
    # gradient/Hessian entry for that dimension -- .nlmSetupEnv() now runs
    # the same .guardScaleC() safety net FOCEi already uses (R/nlmShared.R).
    .fB <- .nlmixr(
      .oneCmtBounded,
      nlmixr2data::theo_sd,
      est = "bobyqa",
      control = bobyqaControl(print = 0L, calcTables = FALSE)
    )
    for (.hm in c("fd", "bfgs", "sr1", "bofill")) {
      .fT <- .nlmixr(
        .oneCmtBounded,
        nlmixr2data::theo_sd,
        est = "trust",
        control = trustControl(print = 0L, calcTables = FALSE, hessianMethod = .hm)
      )
      expect_true(is.finite(.fT$objective), info = .hm)
      expect_equal(.fT$objective, .fB$objective, tolerance = 1e-2, info = .hm)
      expect_equal(unname(.fT$theta), unname(.fB$theta), tolerance = 1e-2, info = .hm)
    }
  })

  test_that("hessianMethod actually changes the Hessian construction", {
    skip_on_cran()
    # A quasi-Newton update (bfgs/sr1/bofill) builds a genuinely different
    # matrix than the fresh finite-difference-of-gradient every call -- this
    # is positive evidence hessianMethod is wired through to nlmTrustFit(),
    # not just that the objective/theta happen to still agree (a well-posed
    # problem can converge to the same point from several different Hessian
    # approximations).
    .rFd <- .nlmixr(
      .oneCmt,
      nlmixr2data::theo_sd,
      est = "trust",
      control = trustControl(print = 0L, calcTables = FALSE, returnTrust = TRUE, hessianMethod = "fd")
    )
    .rSr1 <- .nlmixr(
      .oneCmt,
      nlmixr2data::theo_sd,
      est = "trust",
      control = trustControl(print = 0L, calcTables = FALSE, returnTrust = TRUE, hessianMethod = "sr1")
    )
    expect_true(max(abs(.rFd$hessian - .rSr1$hessian)) > 1)
  })

  test_that("est='trust' returnTrust=TRUE gives the raw trust_solve_c() output", {
    skip_on_cran()
    .raw <- .nlmixr(
      .oneCmt,
      nlmixr2data::theo_sd,
      est = "trust",
      control = trustControl(print = 0L, calcTables = FALSE, returnTrust = TRUE)
    )
    expect_true(is.list(.raw))
    expect_true(is.logical(.raw$converged))
    expect_true(is.numeric(.raw$value))
    expect_equal(dim(.raw$hessian), c(4L, 4L))
  })

  test_that("est='trust' uses its own Hessian only when it is a finite-difference one or asked for (issue 1140)", {
    skip_on_cran()
    .hessCalls <- 0L
    .hess <- nlmixr2Hess
    .fitTrust <- function(...) {
      testthat::with_mocked_bindings(
        nlmixr2Hess = function(...) {
          .hessCalls <<- .hessCalls + 1L
          .hess(...)
        },
        .package = "nlmixr2est",
        .nlmixr(
          .oneCmt,
          nlmixr2data::theo_sd,
          est = "trust",
          control = trustControl(print = 0L, calcTables = FALSE, ...)
        )
      )
    }
    # "trust": the last iteration's (quasi-Newton, by default) Hessian
    .fT <- .fitTrust(covMethod = "trust")
    expect_identical(.hessCalls, 0L)
    expect_true(is.finite(.fT$objective))
    expect_true(.fT$covMethod %in% c("r (trust)", "r+ (trust)", "|r| (trust)"))
    # "r" with hessianMethod = "fd": that Hessian is already a finite-difference one
    .fF <- .fitTrust(hessianMethod = "fd")
    expect_identical(.hessCalls, 0L)
    expect_true(.fF$covMethod %in% c("r", "r+", "|r|"))
    # "r" with the default quasi-Newton hessianMethod: nlmixr2Hess()
    .fR <- .fitTrust()
    expect_identical(.hessCalls, 1L)
    expect_true(.fR$covMethod %in% c("r", "r+", "|r|"))
  })

  test_that("est='trust' warns when the Newton decrement contradicts trust_solve_c()'s own converged flag", {
    # The false-positive-convergence case this warning guards against (RcppTrust
    # reports converged==TRUE off its own internal step-size tolerance while the
    # true Newton step at that point is still large, verified in nlmTrustFit(),
    # src/nlm.cpp) is a genuine numerical edge case, not reliably reproducible
    # on demand from a real fit -- test the warning logic directly instead.
    expect_warning(
      .trustWarnUnderConverged(list(underConverged = TRUE)),
      "Newton step|verified stationary point"
    )
    expect_silent(.trustWarnUnderConverged(list(underConverged = FALSE)))
    expect_silent(.trustWarnUnderConverged(list()))
  })

  test_that("a trust solver error installs no covariance", {
    skip_on_cran()
    # a failed solve hands back a zero-filled Hessian, which used to be
    # regularized by cholSE() into a covariance labelled "r"
    .solve <- nlmTrustFit
    local_mocked_bindings(nlmTrustFit = function(theta, control) {
      .r <- .solve(theta, control)
      .r$error <- -1L
      .r$hessian[] <- 0
      .r
    })
    .fT <- .nlmixr(
      .oneCmt,
      nlmixr2data::theo_sd,
      est = "trust",
      control = trustControl(print = 0L, calcTables = FALSE, covMethod = "trust")
    )
    expect_identical(.fT$covMethod, "failed")
    expect_null(.fT$cov)
    expect_true("R matrix is not finite; covariance step failed" %in% .fT$runInfo)
  })
})
