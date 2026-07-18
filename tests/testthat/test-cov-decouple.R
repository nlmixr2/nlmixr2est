nmTest({
  # Quick, always-run covariance-decoupling checks: the switched family defaults
  # and the setCov(analytic) no-fallback guarantee.  The heavier sa/imp recompute
  # sweep lives in test-cov-decouple-saimp.R (weekly batch).

  .lc <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  .d <- nlmixr2data::theo_sd

  test_that("FOCEI-family default covMethod is now r,s", {
    .f <- suppressWarnings(nlmixr2(.lc, .d, est = "focei",
                                   control = foceiControl(print = 0L)))
    expect_equal(.f$covMethod, "r,s")
    # foceControl / laplaceControl inherit the foceiControl default
    expect_equal(foceiControl()$covMethod, 1L)     # "r,s" -> integer slot 1
    expect_equal(foceControl()$covMethod, 1L)
    expect_equal(laplaceControl()$covMethod, 1L)
  })

  test_that("setCov('analytic') never silently downgrades to r,s", {
    # a linCmt() model is out of analytic-covariance scope
    .f <- suppressWarnings(nlmixr2(.lc, .d, est = "focei",
                                   control = foceiControl(print = 0L)))
    expect_equal(.f$covMethod, "r,s")
    .cov0 <- .f$cov
    # analytic cannot be computed -> error, covariance left unchanged (NOT r,s
    # mislabeled as analytic)
    expect_error(setCov(.f, "analytic"), "could not be computed")
    expect_equal(.f$covMethod, "r,s")
    expect_equal(.f$cov, .cov0)
  })

  test_that("all controls accept sa/imp as a covMethod choice", {
    expect_equal(foceiControl(covMethod = "sa")$covMethodDeferred, "sa")
    expect_equal(foceiControl(covMethod = "imp")$covMethodDeferred, "imp")
    expect_equal(saemControl(covMethod = "imp")$covMethodDeferred, "imp")
    expect_silent(impmapControl(covMethod = "sa"))
    expect_equal(nlmControl(covMethod = "sa")$covMethod, "sa")
    expect_equal(bobyqaControl(covMethod = "imp")$covMethod, "imp")
    # nlme now defaults to keeping its own covariance
    expect_equal(eval(formals(nlmeControl)$covMethod)[1], "nlme")
    # vae now defaults to r,s
    expect_equal(eval(formals(vaeControl)$covMethod)[1], "r,s")
  })
})
