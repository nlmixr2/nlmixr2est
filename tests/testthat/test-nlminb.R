nmTest({
  .pk <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(centr) <- ka * depot - cl / v * centr
      cp <- centr / v
      cp ~ add(add.sd)
    })
  }

  test_that("nlminbControl() defaults to the nlmixr2Hess() covariance, as documented (issue 1140)", {
    expect_identical(nlminbControl()$covMethod, "r")
    expect_identical(nlminbControl(solveType = "grad")$covMethod, "r")
    expect_identical(nlminbControl(covMethod = "nlminb")$covMethod, "nlminb")
    expect_warning(.ctl <- nlminbControl(covMethod = "nlminb", solveType = "fun"), "switching to covMethod='r'")
    expect_identical(.ctl$covMethod, "r")
  })

  test_that("covMethod = \"nlminb\" installs nlminb's own Hessian at the estimates (issue 1140)", {
    skip_on_cran()
    .acc <- new.env(parent = emptyenv())
    .hess <- .nlmixrNlminbHessC
    local_mocked_bindings(
      .nlmixrNlminbHessC = function(pars) {
        .acc$par <- pars + 0
        .acc$h <- .hess(pars)
        .acc$h
      },
      nlmixr2Hess = function(...) stop("nlmixr2Hess() is not nlminb's Hessian")
    )
    .fit <- .nlmixr(
      .pk,
      nlmixr2data::theo_sd,
      est = "nlminb",
      control = nlminbControl(print = 0L, covMethod = "nlminb")
    )
    expect_identical(.fit$covMethod, "r (nlminb)")
    # the last Hessian nlminb's Hessian function gave is at the final (scaled)
    # estimates, and it is the one the covariance is inverted from
    .n <- .fit$env$nlminb
    expect_identical(unname(.acc$par), unname(.n$par.scaled))
    expect_identical(unname(.n$hessian), unname(matrix(.acc$h, 4)))
    expect_equal(unname(.n$cov.scaled), unname(solve(matrix(.acc$h, 4))), tolerance = 1e-8)
  })
})
