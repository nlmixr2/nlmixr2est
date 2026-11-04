# covGrad (#1217): the full R from central differences of the analytic outer gradient

.covGradOne <- function() {
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
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

.covGradFit <- function(...) {
  suppressMessages(nlmixr2(.covGradOne, nlmixr2data::theo_sd, "focei", control = foceiControl(print = 0L, ...)))
}

test_that("covGrad differences the analytic gradient and matches the analytic R", {
  skip_on_cran()
  .g <- .covGradFit(fast = TRUE, covMethod = "r,s")
  expect_true(.g$env$covGrad)
  expect_identical(.g$covMethod, "r,s (full)")
  .np <- nrow(.g$cov)
  expect_identical(.g$env$covEvals[["fullGill"]], 0L)
  # 2p + 1 gradients plus one to four step probes per theta (4 here)
  expect_gte(.g$env$covEvals[["fullR"]], 2L * .np + 1L + 4L)
  expect_lte(.g$env$covEvals[["fullR"]], 2L * .np + 1L + 16L)
  .an <- .covGradFit(fast = TRUE, covMethod = "analytic")
  .r <- setCov(.g, "r (full)")$cov
  .nm <- rownames(.r)
  expect_equal(sqrt(diag(.r)), sqrt(diag(.an$cov[.nm, .nm])), tolerance = 2e-3)
  .o <- .covGradFit(fast = TRUE, covMethod = "r,s", covGrad = FALSE)
  expect_null(.o$env$covGrad)
  expect_gt(.o$env$covEvals[["fullR"]], .g$env$covEvals[["fullR"]])
  expect_equal(sqrt(diag(setCov(.g, "s (full)")$cov)), sqrt(diag(setCov(.o, "s (full)")$cov)), tolerance = 1e-2)
  expect_equal(sqrt(diag(.g$cov)), sqrt(diag(.o$cov)), tolerance = 0.05)
})

test_that("covGrad without the analytic gradient keeps the objective stencil", {
  skip_on_cran()
  .f <- .covGradFit(covMethod = "r,s")
  expect_null(.f$env$covGrad)
  expect_gt(.f$env$covEvals[["fullGill"]], 0L)
})
