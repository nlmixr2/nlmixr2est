nmTest({
  test_that("tidy() p-values are two-sided for positive and negative statistics", {
    skip_if_not_installed("broom.mixed")
    for (fit in list(one.compartment.fit.saem, one.compartment.with.lag.fit.saem)) {
      td <- broom.mixed::tidy(fit, effects = "fixed")
      dfRes <- nobs(fit) - attr(logLik(fit), "df")
      expect_equal(td$p.value, 2 * stats::pt(-abs(td$statistic), dfRes))
      expect_true(all(td$p.value >= 0 & td$p.value <= 1))
    }
    tdLag <- broom.mixed::tidy(one.compartment.with.lag.fit.saem, effects = "fixed")
    expect_true(any(tdLag$statistic < 0))
    expect_true(all(tdLag$p.value[tdLag$statistic < -10] < 1e-10))
  })
})
