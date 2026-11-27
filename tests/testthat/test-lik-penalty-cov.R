# The analytic covariance is the unpenalised observed information, so with a
# likelihood penalty active it declines to the finite-difference Hessian of the
# penalised objective (#1479).

nmTest({
  .lpcMod <- function() {
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
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  .lpcFit <- function(covMethod, ...) {
    suppressMessages(suppressWarnings(nlmixr2(
      .lpcMod,
      nlmixr2data::theo_sd,
      "focei",
      foceiControl(covMethod = covMethod, print = 0L, maxOuterIterations = 0L),
      table = tableControl(cwres = FALSE, npde = FALSE),
      ...
    )))
  }

  test_that("covType=\"analytic\" declines with a likelihood penalty", {
    skip_on_cran()
    .ridge <- list(nlmixr2LikPenaltySpec("thetaL2", thetas = "tcl", par = 50))
    .fa <- .lpcFit("analytic", likPenalty = .ridge)
    .ff <- .lpcFit("r,s", likPenalty = .ridge)
    .f0 <- .lpcFit("analytic")
    expect_identical(.fa$env$covTried$method[1], "analytic")
    expect_identical(.fa$env$covTried$outcome[1], "a likelihood penalty")
    expect_false(identical(.covBaseName(.fa$covMethod), "analytic"))
    expect_equal(.fa$cov, .ff$cov, tolerance = 1e-6)
    # without a penalty the analytic covariance is still used
    expect_true(grepl("analytic", .f0$covMethod, fixed = TRUE))
    expect_false(identical(.f0$env$covTried$method[1], "analytic"))
    expect_equal(.likPenaltyInfo()[["n"]], 0L)
  })
})
