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

  test_that("covMethod = \"n1qn1\" installs n1qn1's quasi-Newton Hessian (issue 1140)", {
    skip_on_cran()
    expect_identical(n1qn1Control()$covMethod, "r")
    local_mocked_bindings(nlmixr2Hess = function(...) stop("nlmixr2Hess() is not n1qn1's Hessian"))
    .fit <- .nlmixr(
      .pk,
      nlmixr2data::theo_sd,
      est = "n1qn1",
      control = n1qn1Control(print = 0L, covMethod = "n1qn1")
    )
    expect_identical(.fit$covMethod, "r (n1qn1)")
    .n <- .fit$env$n1qn1
    expect_identical(unname(.n$hessian), unname(.n$H))
    expect_equal(unname(.n$cov.scaled), unname(solve(.n$H)), tolerance = 1e-8)
  })
})
