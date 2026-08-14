# `lbfgsFactr` tests the objective reduction of a SINGLE step, not stationarity,
# so at the plain 10^-sigdig the analytic-gradient methods stopped as soon as one
# step was small -- 8.6 OFV units short on a 2-cmt oral fit, in 6 outer
# evaluations.  The default is two orders tighter for that reason
# (plans/foceif-outer-opt-pgtol.md).  Control-level checks only, so this stays in
# the push/PR subset; the fit-based check is in test-focei-factr-fit.R.

nmTest({
  test_that("lbfgsFactr is two orders tighter than the other sigdig tolerances", {
    .ctl <- foceiControl()
    expect_equal(.ctl$lbfgsFactr,
                 10^(-.ctl$sigdig - 2) / .Machine$double.eps)
    # the other sigdig-derived optimizer tolerances are NOT shifted -- nlminb's
    # in particular, where tightening measured no benefit
    expect_equal(.ctl$rel.tol, 10^(-.ctl$sigdig))
    expect_equal(.ctl$x.tol, 10^(-.ctl$sigdig))
    expect_equal(.ctl$rhoend, 10^(-.ctl$sigdig))

    # and it tracks sigdig
    expect_equal(foceiControl(sigdig = 5)$lbfgsFactr,
                 1e-7 / .Machine$double.eps)

    # an explicit value still wins
    expect_equal(foceiControl(lbfgsFactr = 1e7)$lbfgsFactr, 1e7)
  })
})
