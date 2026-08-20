nmTest({
  test_that("foceiControl(innerOpt=) trust mapping", {
    expect_equal(foceiControl()$innerOpt, 1L)
    expect_equal(foceiControl(innerOpt = "n1qn1")$innerOpt, 1L)
    expect_equal(foceiControl(innerOpt = "BFGS")$innerOpt, 2L)
    expect_equal(foceiControl(innerOpt = "trust")$innerOpt, 3L)
    expect_equal(foceiControl(innerOpt = 3L)$innerOpt, 3L)
    expect_error(foceiControl(innerOpt = "bogus"))

    expect_equal(foceiControl()$trustConf, 0.975)
    expect_equal(foceiControl(trustConf = 0.9)$trustConf, 0.9)
    expect_null(foceiControl()$trustRinit)
    expect_null(foceiControl()$trustRmax)
    expect_equal(foceiControl(trustRmax = 2)$trustRmax, 2)
    expect_error(foceiControl(trustRinit = 5, trustRmax = 1))

    .ctl <- foceiControl(innerOpt = "trust", trustConf = 0.9)
    expect_equal(do.call(foceiControl, .ctl)$innerOpt, 3L)
    expect_equal(do.call(foceiControl, .ctl)$trustConf, 0.9)
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

  test_that("innerOpt='trust' restart cascade completes when the inner solve doesn't converge", {
    skip_on_cran()
    # A tiny maxInnerIterations forces trust_solve_c() to hit iterlim before
    # converging, exercising the !converged nudge-restart cascade in the
    # trustInner branch of innerOpt1() (src/inner.cpp). Regression coverage
    # for a bug an outside review caught: trustSolveAt() must reset
    # fInd->badSolve per attempt (n1qn1's cascade already does this on every
    # restart) -- without that reset, an NA during the first attempt
    # permanently poisons every later nudge via trustInnerObjfun's own
    # badSolve guard, silently turning the whole cascade into a no-op.
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, est = "focei",
              control = foceiControl(innerOpt = "trust", maxOuterIterations = 5,
                                      maxInnerIterations = 2,
                                      covMethod = "", calcTables = FALSE, print = 0))
    ))
    expect_true(is.finite(.fit$objf))
    expect_true(.nTrustInner() > 0L)
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
