## Population likelihood penalties that depend on omega (#1481): the built-in
## omegaL2 penalty is lambda * sum((omega[k,k] - target)^2) on the -2LL scale.

.likPenOmegaOde <- function() {
  ini({
    tka <- 0.45; tcl <- 1; tv <- 3.45
    eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

test_that("omegaL2 value and the finite-difference omega gradient", {
  .om <- matrix(c(0.4, 0.1, 0.1, 0.3), 2, 2)
  .ptr <- .nlmixr2LikPenaltyBuiltin("omegaL2")
  .Call(
    `_nlmixr2est_likPenaltyActivate`,
    list(.ptr),
    list(integer(0)),
    list(c(10, 0.05, 1)),
    list(NULL),
    NA_integer_
  )
  on.exit(.likPenaltyDeactivate())
  expect_equal(.Call(`_nlmixr2est_likPenaltyPop`, c(1, 2), .om, 10L), 10 * 0.35^2)
  expect_equal(.Call(`_nlmixr2est_likPenaltyPop`, c(1, 2), NULL, 10L), 0)
  .g <- .Call(`_nlmixr2est_likPenaltyPopOmegaGrad`, c(1, 2), .om, 10L)
  expect_equal(.g, matrix(c(20 * 0.35, 0, 0, 0), 2, 2), tolerance = 1e-6)
  # every eta when none is selected
  .Call(
    `_nlmixr2est_likPenaltyActivate`,
    list(.ptr),
    list(integer(0)),
    list(c(1, 0)),
    list(NULL),
    NA_integer_
  )
  expect_equal(.Call(`_nlmixr2est_likPenaltyPop`, 1, .om, 10L), 0.4^2 + 0.3^2)
  .g <- .Call(`_nlmixr2est_likPenaltyPopOmegaGrad`, 1, .om, 10L)
  expect_equal(.g, diag(c(0.8, 0.6)), tolerance = 1e-6)
})

test_that("the focei analytic gradient includes the omega penalty", {
  skip_on_cran()
  .ctl <- foceiControl(print = 0L, covMethod = "", fast = TRUE)
  .f0 <- suppressWarnings(nlmixr2(.likPenOmegaOde, nlmixr2data::theo_sd, est = "focei", control = .ctl))
  .mp <- nlmixr2LikPenaltyAdd(
    .likPenOmegaOde,
    nlmixr2LikPenaltySpec("omegaL2", par = c(10, 0.05, 1), inObjf = TRUE)
  )
  .f1 <- suppressWarnings(nlmixr2(.mp, nlmixr2data::theo_sd, est = "focei", control = .ctl))
  expect_true(.f1$env$nAnalyticGradDirect > 0)
  # the penalty pulls omega[1,1] toward its target
  expect_lt(.f1$omega[1, 1], .f0$omega[1, 1] - 0.05)
  .pen <- 10 * (.f1$omega[1, 1] - 0.05)^2
  expect_equal(.f1$likPenalty$population, .pen, tolerance = 1e-6)
  # at a common point the gradients differ by exactly d(penalty)/d(omega)
  .g0 <- nlmixr2ObjGrad(.f0)
  .g1 <- nlmixr2ObjGrad(.f1, theta = .f0$theta, omega = .f0$omega)
  expect_equal(.g1$objective - .g0$objective, 10 * (.f0$omega[1, 1] - 0.05)^2, tolerance = 1e-3)
  .d <- .g1$gradient - .g0$gradient
  expect_equal(.d[["eta.ka"]], 20 * (.f0$omega[1, 1] - 0.05), tolerance = 0.01)
  expect_lt(max(abs(.d[c("eta.cl", "eta.v")])), 0.05)
})

test_that("saem's omega M-step and objective see the omega penalty", {
  skip_on_cran()
  .pen <- list(nlmixr2LikPenaltySpec("omegaL2", par = c(1000, 0.05, 1)))
  .fit <- function(est, ...) {
    suppressWarnings(nlmixr2(.likPenOmegaOde, nlmixr2data::theo_sd, est = est, control = list(print = 0L), ...))
  }
  .s0 <- .fit("saem")
  .s1 <- .fit("saem", likPenalty = .pen)
  .f1 <- .fit("focei", likPenalty = .pen)
  # the penalty pulls omega[1,1] toward its target under both methods, and
  # they agree within Monte Carlo error
  expect_lt(.s1$omega[1, 1], .s0$omega[1, 1] - 0.2)
  expect_lt(abs(.s1$omega[1, 1] - .f1$omega[1, 1]), 0.05)
  # the stored value is the one in the objective: inObjf = FALSE takes back
  # exactly the penalty the objective held
  expect_equal(.s1$likPenalty$population, 1000 * (.s1$omega[1, 1] - 0.05)^2, tolerance = 1e-6)
  .s2 <- .fit("saem", likPenalty = list(nlmixr2LikPenaltySpec("omegaL2", par = c(1000, 0.05, 1), inObjf = TRUE)))
  expect_equal(.s2$objf - .s1$objf, .s1$likPenalty$population, tolerance = 1e-4)
})

test_that("a theta-only penalty leaves saem's omega update unchanged", {
  skip_on_cran()
  .ctl <- saemControl(print = 0L, seed = 42)
  .p0 <- list(nlmixr2LikPenaltySpec("thetaL2", thetas = "tka", par = c(0, 0, 100), inObjf = TRUE))
  .s0 <- suppressWarnings(nlmixr2(.likPenOmegaOde, nlmixr2data::theo_sd, est = "saem", control = .ctl))
  .s1 <- suppressWarnings(nlmixr2(
    .likPenOmegaOde,
    nlmixr2data::theo_sd,
    est = "saem",
    control = .ctl,
    likPenalty = .p0
  ))
  expect_equal(.s1$omega, .s0$omega)
  expect_equal(.s1$objf - .s0$objf, 100, tolerance = 1e-6)
})
