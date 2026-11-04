# sigdig sets the ODE rtol to 10^(-2*sigdig) (issue #1131): a 10^-sigdig optimizer
# needs an objective accurate to its square, and the former rtol = 10^-sigdig
# biased the FOCEi objective by several units.  Control-level only.

nmTest({
  .tols <- function(rx) {
    vapply(
      c("rtol", "atol", "rtolSens", "atolSens", "ssRtol", "ssAtol"),
      function(f) max(rx[[f]]), numeric(1)
    )
  }
  .expected <- function(rtol) {
    c(
      rtol = rtol, atol = rtol / 100, rtolSens = rtol, atolSens = rtol / 100,
      ssRtol = 10 * rtol, ssAtol = rtol / 10
    )
  }

  test_that("default foceiControl solves at rxode2's own 1e-6/1e-8", {
    expect_equal(.tols(foceiControl()$rxControl), .expected(1e-6))
    expect_equal(.tols(foceiControl(sigdig = 4)$rxControl), .expected(1e-8))
    expect_equal(
      foceiControl()$newuoaRhoend,
      10^(-foceiControl()$sigdig)
    )
  })

  test_that("the sigdig ODE rtol is floored at 1e-12", {
    expect_equal(.tols(foceiControl(sigdig = 8)$rxControl), .expected(1e-12))
  })

  test_that("other sigdig-scaled methods share the mapping", {
    expect_equal(.tols(bobyqaControl()$rxControl), .expected(1e-6))
    expect_equal(.tols(nlmControl()$rxControl), .expected(1e-6))
    expect_equal(.tols(nlsControl()$rxControl), .expected(1e-8))
  })

  test_that("an explicit tolerance in a rxControl list still wins", {
    .rx <- foceiControl(rxControl = list(rtol = 1e-4))$rxControl
    expect_equal(max(.rx$rtol), 1e-4)
    expect_equal(max(.rx$atol), 1e-8)
  })
})
