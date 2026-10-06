nmTest({
  test_that("foceiControl() innerLbfgs* controls", {
    expect_equal(foceiControl(innerOpt = "BFGS")$innerOpt, 2L)
    expect_equal(foceiControl()$innerLbfgsLmm, 5L)
    expect_equal(foceiControl(innerLbfgsLmm = 3)$innerLbfgsLmm, 3L)
    expect_error(foceiControl(innerLbfgsLmm = 0))
    expect_error(foceiControl(innerLbfgsLmm = 2.5))
    expect_equal(foceiControl(sigdig = 4)$innerLbfgsFactr, 1e-6 / .Machine$double.eps)
    expect_equal(foceiControl(sigdig = 14)$innerLbfgsFactr, 1)
    expect_equal(foceiControl()$innerLbfgsPgtol, 0)
    expect_equal(foceiControl(sigdig = 3)$innerLbfgsAbstol, 1e-5)
    expect_equal(foceiControl(sigdig = 3)$innerLbfgsReltol, 1e-5)
    # separate from the outer L-BFGS-B tolerances
    expect_equal(foceiControl(sigdig = 3, abstol = 0.1)$innerLbfgsAbstol, 1e-5)
    expect_equal(foceiControl(innerLbfgsFactr = 10)$lbfgsFactr, foceiControl()$lbfgsFactr)
    expect_error(foceiControl(innerLbfgsPgtol = -1))
    expect_error(foceiControl(innerLbfgsAbstol = -1))
    expect_error(foceiControl(innerLbfgsReltol = Inf))
    .ctl <- foceiControl(innerOpt = "BFGS", innerLbfgsLmm = 4L, innerLbfgsFactr = 100)
    expect_equal(do.call(foceiControl, .ctl)$innerLbfgsLmm, 4L)
    expect_equal(do.call(foceiControl, .ctl)$innerLbfgsFactr, 100)
  })

  test_that("innerOpt='BFGS' errors clearly with an old lbfgsb3c", {
    expect_error(.foceiAssertInnerBfgs(2L, have = FALSE), "lbfgsb3c >= 2024-3.6")
    expect_true(.foceiAssertInnerBfgs(2L, have = TRUE))
    expect_true(.foceiAssertInnerBfgs(1L, have = FALSE))
    expect_true(.lbfgsb3ctsAvailable())
  })

  .oneCmtBfgs <- function() {
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

  .fitBfgsCmp <- function(innerOpt, cores = 1L, ...) {
    suppressWarnings(suppressMessages(
      nlmixr2(
        .oneCmtBfgs,
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

  test_that("innerOpt='BFGS' runs L-BFGS-B and matches n1qn1/trust", {
    skip_on_cran()
    .fb <- .fitBfgsCmp("BFGS")
    .fn <- .fitBfgsCmp("n1qn1")
    .ft <- .fitBfgsCmp("trust")
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

  test_that("innerOpt='BFGS' gives identical fits with cores=1 and cores=2", {
    skip_on_cran()
    .f1 <- .fitBfgsCmp("BFGS", cores = 1L)
    .f2 <- .fitBfgsCmp("BFGS", cores = 2L)
    expect_gt(.f2$env$nLbfgsInner[["calls"]], 0L)
    expect_identical(.f1$objf, .f2$objf)
    expect_identical(as.data.frame(.f1$eta), as.data.frame(.f2$eta))
    expect_identical(.f1$theta, .f2$theta)
  })

  test_that("innerOpt='BFGS' with mceta restarts and a dnorm() endpoint", {
    skip_on_cran()
    .fm <- .fitBfgsCmp("BFGS", mceta = 3L)
    expect_true(is.finite(.fm$objf))
    expect_gt(.fm$env$nLbfgsInner[["calls"]], 0L)
    .ll <- .oneCmtBfgs |> model(linCmt() ~ add(add.sd) + dnorm())
    .fl <- suppressWarnings(suppressMessages(
      nlmixr2(.ll, nlmixr2data::theo_sd, est = "focei",
        control = foceiControl(innerOpt = "BFGS", covMethod = "",
          calcTables = FALSE, print = 0))
    ))
    .fn <- suppressWarnings(suppressMessages(
      nlmixr2(.ll, nlmixr2data::theo_sd, est = "focei",
        control = foceiControl(innerOpt = "n1qn1", covMethod = "",
          calcTables = FALSE, print = 0))
    ))
    expect_gt(.fl$env$nLbfgsInner[["calls"]], 0L)
    expect_equal(.fl$objf, .fn$objf, tolerance = 1e-3)
  })

  test_that("innerOpt='BFGS' works for est='laplace'", {
    skip_on_cran()
    .fl <- suppressWarnings(suppressMessages(
      nlmixr2(.oneCmtBfgs, nlmixr2data::theo_sd, est = "laplace",
        control = laplaceControl(innerOpt = "BFGS", covMethod = "",
          calcTables = FALSE, print = 0))
    ))
    expect_true(is.finite(.fl$objf))
    expect_gt(.fl$env$nLbfgsInner[["calls"]], 0L)
  })

  test_that("innerOpt='BFGS' counts maxit stops", {
    skip_on_cran()
    .f <- .fitBfgsCmp("BFGS", maxInnerIterations = 2L, maxOuterIterations = 2L)
    .n <- .f$env$nLbfgsInner
    expect_gt(.n[["maxit"]], 0L)
    expect_lte(.n[["maxit"]], .n[["notConverged"]])
  })
})
