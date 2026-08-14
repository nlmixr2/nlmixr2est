# L-BFGS-B's `factr` tests the objective reduction of a SINGLE step, so with
# pgtol = 0 it stops when one step is small rather than when the gradient is
# flat.  `pgtol` is the stationarity test that catches it, but only an analytic
# gradient can support it.  These are control-level checks only, so this stays in
# the push/PR subset; the fit-based convergence check lives in
# test-focei-pgtol-fit.R (weekly-batched).

nmTest({
  test_that("the outer projected-gradient test follows the analytic gradient", {
    # fast=FALSE -> finite-difference outer gradient -> suppressed
    .ctl <- foceiControl()
    expect_equal(.ctl$lbfgsPgtolOuter, 0)
    expect_true(.ctl$lbfgsPgtolOuterDefault)

    # fast=TRUE -> analytic outer gradient -> derived from sigdig
    .ctl <- foceiControl(fast = TRUE)
    expect_equal(.ctl$lbfgsPgtolOuter, 10^(-.ctl$sigdig))
    .ctl <- foceiControl(fast = TRUE, sigdig = 6)
    expect_equal(.ctl$lbfgsPgtolOuter, 1e-6)

    # the *f wrappers get it through fast=TRUE
    .ctl <- getValidNlmixrCtl(structure(list(NULL),
                                        class = c("foceif", "getValidNlmixrControl")))
    expect_true(.ctl$lbfgsPgtolOuter > 0)

    # a derivative-free outerOpt reverts fast, which must re-suppress it
    expect_warning(.ctl <- foceiControl(fast = TRUE, outerOpt = "bobyqa"),
                   "derivative-free")
    expect_equal(.ctl$lbfgsPgtolOuter, 0)
  })

  test_that("the inner projected-gradient test is always derived", {
    # the inner eta gradient is always the sensitivity gradient, so unlike the
    # outer one it does not depend on `fast`
    expect_equal(foceiControl()$lbfgsPgtol, 10^(-foceiControl()$sigdig))
    expect_equal(foceiControl(fast = TRUE)$lbfgsPgtol, 10^(-foceiControl()$sigdig))
    expect_equal(foceiControl(sigdig = 5)$lbfgsPgtol, 1e-5)
  })

  test_that("an explicit lbfgsPgtolOuter survives a control round-trip", {
    .ctl <- foceiControl(lbfgsPgtolOuter = 1e-9)
    expect_equal(.ctl$lbfgsPgtolOuter, 1e-9)
    expect_false(.ctl$lbfgsPgtolOuterDefault)
    # the *f wrapper rebuilds the control under fast=TRUE; a user value must not
    # be re-derived away, while a defaulted one must be
    .f <- getValidNlmixrCtl(structure(list(.ctl),
                                      class = c("foceif", "getValidNlmixrControl")))
    expect_equal(.f$lbfgsPgtolOuter, 1e-9)
    .f <- getValidNlmixrCtl(structure(list(foceiControl()),
                                      class = c("foceif", "getValidNlmixrControl")))
    expect_true(.f$lbfgsPgtolOuter > 0)
  })

  test_that("the nlm family derives pgtol only where it passes an analytic gradient", {
    # est="lbfgsb3c" always passes gr=.nlmixrOptimGradC
    expect_equal(lbfgsb3cControl()$pgtol, 10^(-lbfgsb3cControl()$sigdig))
    expect_equal(lbfgsb3cControl(sigdig = 5)$pgtol, 1e-5)

    # est="optim" only passes gr for method="L-BFGS-B" with solveType="grad"
    expect_equal(optimControl(method = "L-BFGS-B")$pgtol,
                 10^(-optimControl()$sigdig))
    expect_equal(optimControl(method = "L-BFGS-B", solveType = "fun")$pgtol, 0)
    expect_equal(optimControl(method = "BFGS")$pgtol, 0)
    expect_equal(optimControl()$pgtol, 0) # Nelder-Mead
  })
})
