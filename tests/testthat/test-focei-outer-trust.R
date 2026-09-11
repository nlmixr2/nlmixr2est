test_that("outerOpt='trust' validates its own control", {
  expect_error(foceiControl(outerOpt = "trust", outerTrustHessian = "analytic"),
               "requires fast")
  expect_no_error(foceiControl(outerOpt = "trust", fast = TRUE,
                               outerTrustHessian = "analytic"))
  expect_error(foceiControl(outerTrustRinit = 0), "must be > 0")
  expect_error(foceiControl(outerTrustRmax = 0), "must be > 0")
  expect_error(foceiControl(outerTrustRinit = 2, outerTrustRmax = 1),
               "cannot be larger")
  expect_error(foceiControl(outerTrustRelStep = 0), "must be > 0")
  expect_error(foceiControl(outerTrustRestarts = -1L))
  # a derivative-based outer optimizer must NOT trip the derivative-free
  # fast= downgrade that bobyqa/uobyqa/newuoa do
  expect_true(foceiControl(outerOpt = "trust", fast = TRUE)$fast)
  expect_equal(foceiControl(outerOpt = "trust")$outerOptTxt, "trust")
})

test_that("the Newton decrement gate reads a trust result", {
  # positive definite Hessian: 0.5 * g' H^-1 g
  expect_equal(.trustOuterDecrement(list(gradient = c(1, 2),
                                         hessian = diag(c(2, 8)))),
               0.5 * (1 / 2 + 4 / 8))
  # an indefinite Hessian at a reported minimum is not a minimum
  expect_true(is.na(.trustOuterDecrement(list(gradient = c(1, 0),
                                              hessian = diag(c(-1, 1))))))
  expect_true(is.na(.trustOuterDecrement(list(gradient = c(NA_real_, 0),
                                              hessian = diag(2)))))
  expect_true(is.na(.trustOuterDecrement(list(gradient = c(1, 0),
                                              hessian = NULL))))
})

test_that("outerOpt='trust' fits and consumes the analytic outer Hessian", {
  skip_on_cran()
  model <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
          eta.cl ~ 0.3; add.sd <- 0.7 })
    model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
            d/dt(depot) <- -ka * depot
            d/dt(center) <- ka * depot - cl / v * center
            cp <- center / v
            cp ~ add(add.sd) })
  }
  d <- nlmixr2data::theo_sd
  ctl <- function(...) {
    foceiControl(print = 0L, calcTables = FALSE, covMethod = "",
                 outerOpt = "trust", ...)
  }
  fitA <- .nlmixr(model, d, "focei", ctl(fast = TRUE))
  # The counter is what proves the analytic Hessian ran; equal objectives alone
  # cannot tell it from the quasi-Newton fallback.
  expect_gt(fitA$env$optReturn$hessianEvaluations, 0L)
  expect_false(fitA$env$optReturn$hessianFallback)
  expect_true(is.finite(fitA$objf))

  fitB <- .nlmixr(model, d, "focei", ctl(fast = TRUE, outerTrustHessian = "bfgs"))
  expect_equal(fitB$env$optReturn$hessianEvaluations, 0L)
  expect_equal(fitB$objf, fitA$objf, tolerance = 1e-3)

  # ... and it reaches the same optimum as the shipping outer optimizer
  fitN <- .nlmixr(model, d, "focei",
                  foceiControl(print = 0L, calcTables = FALSE, covMethod = "",
                               outerOpt = "nlminb", fast = TRUE))
  expect_equal(fitA$objf, fitN$objf, tolerance = 1e-3)
  expect_equal(unname(fixef(fitA)), unname(fixef(fitN)), tolerance = 1e-2)
})
