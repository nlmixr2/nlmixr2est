# Fit-based check of the same()-block registry (.iovSameMethods and
# .sameBlockMethods, R/iov.R): every method it lets through really keeps a
# user's repeated block.  The refusal of the others is in test-iov-same.R.
# Slow (one fit per method); weekly via .slowBatches in tests/testthat.R.

.sameFitMod <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ same()
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    # an ODE, since the full conditional-Hessian methods refuse linCmt()
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

.sameFitOmega <- function(fit) {
  .om <- fit$omega
  if (is.list(.om)) .om$id else .om
}

test_that("every estimating method the registry lets through keeps eta.v tied to eta.cl", {
  skip_on_cran()
  .ests <- c(.iovSameMethods, names(.sameBlockMethods)[.sameBlockMethods == "keeps the tie, measured"])
  # the 33 FOCEi-family methods and the five measured to keep the tie
  expect_length(.ests, 38L)
  for (.e in .ests) {
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(.sameFitMod, nlmixr2data::theo_sd, .e, control = list(print = 0))
    ))
    .om <- .sameFitOmega(.fit)
    expect_equal(.om["eta.v", "eta.v"], .om["eta.cl", "eta.cl"], tolerance = 1e-12, label = .e)
    # the variance was estimated, so the tie is not just the shared initial value
    expect_false(isTRUE(all.equal(.om["eta.cl", "eta.cl"], 0.3)), label = .e)
  }
})

test_that("the methods that do not estimate omega run on a same() model", {
  skip_on_cran()
  .ests <- names(.sameBlockMethods)[.sameBlockMethods == "does not estimate omega"]
  # "output" builds a fit's output from an estimate; nlmixr2() does not dispatch to it
  .ests <- setdiff(.ests, "output")
  expect_identical(.ests, c("posthoc", "predict", "rxSolve", "simulate", "simulation"))
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(.sameFitMod, nlmixr2data::theo_sd, "posthoc", control = list(print = 0))
  ))
  .om <- .sameFitOmega(.fit)
  expect_equal(unname(diag(.om)[c("eta.cl", "eta.v")]), c(0.3, 0.3))
  for (.e in setdiff(.ests, "posthoc")) {
    .sim <- suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(.sameFitMod, nlmixr2data::theo_sd, .e)
    ))
    expect_s3_class(.sim, "rxSolve")
    expect_gt(nrow(.sim), 0L, label = .e)
  }
})

test_that("foceiLikLoad() loads a same() model for an est with no method", {
  skip_on_cran()
  # the hooks run with the est foceiLikLoad() is given; with no nlmixr2Est
  # method nothing refuses the block, and the FOCEi inner problem keeps the tie
  .h <- suppressMessages(foceiLikLoad(.sameFitMod, nlmixr2data::theo_sd, "focei", est = "x"))
  withr::defer(foceiLikUnload())
  expect_identical(.h$neta, 3L)
})
