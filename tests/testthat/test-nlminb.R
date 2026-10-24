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

  test_that("nlminbControl() defaults to nlminb's own Hessian only where the fit computes it (issue 1140)", {
    expect_identical(nlminbControl()$covMethod, "nlminb")
    expect_identical(nlminbControl(solveType = "hessian")$covMethod, "nlminb")
    expect_identical(nlminbControl(solveType = "grad")$covMethod, "r")
    expect_identical(nlminbControl(solveType = "fun")$covMethod, "r")
    expect_identical(nlminbControl(covMethod = "r")$covMethod, "r")
    expect_identical(nlminbControl(covMethod = "nlminb")$covMethod, "nlminb")
    expect_warning(.ctl <- nlminbControl(covMethod = "nlminb", solveType = "fun"), "switching to covMethod='r'")
    expect_identical(.ctl$covMethod, "r")
  })

  test_that("the nlm-family covariance Hessian runs at the covariance probe tolerances (issue 1140)", {
    skip_on_cran()
    # at the fit's own rtol (1e-3) a "fun" fit's stencil differenced solver noise:
    # an indefinite R, repaired to |r| with SEs about 6x too small
    .fun <- .nlmixr(.pk, nlmixr2data::theo_sd, est = "nlminb", control = nlminbControl(print = 0L, solveType = "fun"))
    .grad <- .nlmixr(.pk, nlmixr2data::theo_sd, est = "nlminb", control = nlminbControl(print = 0L, solveType = "grad"))
    expect_identical(.fun$covMethod, "r")
    expect_identical(.grad$covMethod, "r")
    expect_equal(sqrt(diag(.fun$cov)), sqrt(diag(.grad$cov)), tolerance = 0.05)
  })

  test_that("covMethod = \"nlminb\" differences the analytic gradient at the estimates (issue 1140)", {
    skip_on_cran()
    .acc <- new.env(parent = emptyenv())
    .gh <- .nlmGradHessian
    local_mocked_bindings(
      .nlmGradHessian = function(par) {
        .acc$par <- par + 0
        .acc$h <- .gh(par)
        # an independent central difference of the analytic gradient there
        .acc$fd <- vapply(
          seq_along(.acc$par),
          function(k) {
            .e <- replace(numeric(length(.acc$par)), k, 1e-4)
            (.nlmixrNlminbGradC(.acc$par + .e) - .nlmixrNlminbGradC(.acc$par - .e)) / 2e-4
          },
          numeric(length(.acc$par))
        )
        .acc$h
      },
      nlmixr2Hess = function(...) stop("nlmixr2Hess() is not the \"nlminb\" Hessian")
    )
    for (.st in c("hessian", "grad")) {
      .acc$par <- NULL
      .fit <- .nlmixr(
        .pk,
        nlmixr2data::theo_sd,
        est = "nlminb",
        control = nlminbControl(print = 0L, covMethod = "nlminb", solveType = .st)
      )
      expect_identical(.fit$covMethod, "r (nlminb)", info = .st)
      .n <- .fit$env$nlminb
      expect_identical(unname(.acc$par), unname(.n$par.scaled), info = .st)
      expect_identical(unname(.n$hessian), unname(.acc$h), info = .st)
      expect_equal(unname(.n$cov.scaled), unname(solve(.acc$h)), tolerance = 1e-8, info = .st)
      expect_equal(.acc$h, 0.5 * (.acc$fd + t(.acc$fd)), tolerance = 1e-3, info = .st)
    }
  })
})
