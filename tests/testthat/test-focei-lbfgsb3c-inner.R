nmTest({
  test_that("foceiControl() innerLbfgs* controls", {
    expect_equal(foceiControl(innerOpt = "lbfgsb3c")$innerOpt, 2L)
    expect_equal(foceiControl()$innerLbfgsLmm, 5L)
    expect_equal(foceiControl(innerLbfgsLmm = 3)$innerLbfgsLmm, 3L)
    expect_error(foceiControl(innerLbfgsLmm = 0))
    expect_error(foceiControl(innerLbfgsLmm = 2.5))
    expect_equal(foceiControl(sigdig = 4)$innerLbfgsFactr, 1e-6 / .Machine$double.eps)
    expect_equal(foceiControl(sigdig = 14)$innerLbfgsFactr, 1)
    # every inner tolerance follows sigdig
    expect_equal(foceiControl(sigdig = 3)$innerLbfgsPgtol, 1e-5)
    expect_equal(foceiControl(sigdig = 5)$innerLbfgsPgtol, 1e-7)
    expect_equal(foceiControl(sigdig = 5)$innerLbfgsAbstol, 1e-7)
    expect_equal(foceiControl(sigdig = 5)$innerLbfgsFactr, 1e-7 / .Machine$double.eps)
    expect_equal(foceiControl(innerLbfgsPgtol = 0)$innerLbfgsPgtol, 0)
    expect_equal(foceiControl(sigdig = 3)$innerLbfgsAbstol, 1e-5)
    expect_equal(foceiControl(sigdig = 3)$innerLbfgsReltol, 1e-5)
    # separate from the outer L-BFGS-B tolerances
    expect_equal(foceiControl(sigdig = 3, abstol = 0.1)$innerLbfgsAbstol, 1e-5)
    expect_equal(foceiControl(innerLbfgsFactr = 10)$lbfgsFactr, foceiControl()$lbfgsFactr)
    expect_error(foceiControl(innerLbfgsPgtol = -1))
    expect_error(foceiControl(innerLbfgsAbstol = -1))
    expect_error(foceiControl(innerLbfgsReltol = Inf))
    .ctl <- foceiControl(innerOpt = "lbfgsb3c", innerLbfgsLmm = 4L, innerLbfgsFactr = 100)
    expect_equal(do.call(foceiControl, .ctl)$innerLbfgsLmm, 4L)
    expect_equal(do.call(foceiControl, .ctl)$innerLbfgsFactr, 100)
  })

  test_that("innerOpt='lbfgsb3c' errors clearly with an old lbfgsb3c", {
    expect_error(.foceiAssertInnerLbfgsb3c(2L, have = FALSE), "lbfgsb3c >= 2024-3.6")
    expect_true(.foceiAssertInnerLbfgsb3c(2L, have = TRUE))
    expect_true(.foceiAssertInnerLbfgsb3c(1L, have = FALSE))
    expect_true(.lbfgsb3ctsAvailable())
  })

  .oneCmtLbfgsb3c <- function() {
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
      linCmt() ~ add(add.sd)
    })
  }

  .fitLbfgsb3cCmp <- function(innerOpt, cores = 1L, ...) {
    suppressWarnings(suppressMessages(
      nlmixr2(
        .oneCmtLbfgsb3c,
        nlmixr2data::theo_sd,
        est = "focei",
        control = foceiControl(
          innerOpt = innerOpt,
          covMethod = "",
          calcTables = FALSE,
          print = 0,
          rxControl = rxode2::rxControl(cores = cores),
          ...
        )
      )
    ))
  }

  test_that("innerOpt='lbfgsb3c' runs L-BFGS-B and matches n1qn1/trust", {
    skip_on_cran()
    .fb <- .fitLbfgsb3cCmp("lbfgsb3c")
    .fn <- .fitLbfgsb3cCmp("n1qn1")
    .ft <- .fitLbfgsb3cCmp("trust")
    expect_true(is.finite(.fb$objf))
    expect_equal(.fb$objf, .fn$objf, tolerance = 1e-3)
    expect_equal(.fb$objf, .ft$objf, tolerance = 1e-3)
    expect_equal(as.data.frame(.fb$eta), as.data.frame(.fn$eta), tolerance = 5e-2)
    # The L-BFGS-B path actually ran, and converged.
    .n <- .fb$env$nLbfgsInner
    expect_gt(.n[["calls"]], 0L)
    expect_equal(.n[["notConverged"]], 0L)
    expect_null(.fn$env$nLbfgsInner)
  })

  test_that("innerOpt='lbfgsb3c' gives identical fits with cores=1 and cores=2", {
    skip_on_cran()
    .f1 <- .fitLbfgsb3cCmp("lbfgsb3c", cores = 1L)
    .f2 <- .fitLbfgsb3cCmp("lbfgsb3c", cores = 2L)
    expect_gt(.f2$env$nLbfgsInner[["calls"]], 0L)
    expect_identical(.f1$objf, .f2$objf)
    expect_identical(as.data.frame(.f1$eta), as.data.frame(.f2$eta))
    expect_identical(.f1$theta, .f2$theta)
  })

  test_that("innerOpt='lbfgsb3c' with mceta restarts and a dnorm() endpoint", {
    skip_on_cran()
    .fm <- .fitLbfgsb3cCmp("lbfgsb3c", mceta = 3L)
    expect_true(is.finite(.fm$objf))
    expect_gt(.fm$env$nLbfgsInner[["calls"]], 0L)
    .ll <- .oneCmtLbfgsb3c |> model(linCmt() ~ add(add.sd) + dnorm())
    .fl <- suppressWarnings(suppressMessages(
      nlmixr2(
        .ll,
        nlmixr2data::theo_sd,
        est = "focei",
        control = foceiControl(innerOpt = "lbfgsb3c", covMethod = "", calcTables = FALSE, print = 0)
      )
    ))
    .fn <- suppressWarnings(suppressMessages(
      nlmixr2(
        .ll,
        nlmixr2data::theo_sd,
        est = "focei",
        control = foceiControl(innerOpt = "n1qn1", covMethod = "", calcTables = FALSE, print = 0)
      )
    ))
    expect_gt(.fl$env$nLbfgsInner[["calls"]], 0L)
    expect_equal(.fl$objf, .fn$objf, tolerance = 1e-3)
  })

  test_that("innerOpt='lbfgsb3c' works for est='laplace'", {
    skip_on_cran()
    .fl <- suppressWarnings(suppressMessages(
      nlmixr2(
        .oneCmtLbfgsb3c,
        nlmixr2data::theo_sd,
        est = "laplace",
        control = laplaceControl(innerOpt = "lbfgsb3c", covMethod = "", calcTables = FALSE, print = 0)
      )
    ))
    expect_true(is.finite(.fl$objf))
    expect_gt(.fl$env$nLbfgsInner[["calls"]], 0L)
  })

  test_that("innerOpt='lbfgsb3c' counts maxit stops", {
    skip_on_cran()
    .f <- .fitLbfgsb3cCmp("lbfgsb3c", maxInnerIterations = 2L, maxOuterIterations = 2L)
    .n <- .f$env$nLbfgsInner
    expect_gt(.n[["maxit"]], 0L)
    expect_lte(.n[["maxit"]], .n[["notConverged"]])
  })
})
