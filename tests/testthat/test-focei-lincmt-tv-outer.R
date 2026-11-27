nmTest({
  # A time-varying covariate on a linCmt() parameter: the outer gradient's columns come
  # from rxode2lincmt's direction carry (.foceiOuterLinCmtCarry) instead of the
  # parameter-space reads, which read NA there and sent every subject to the
  # per-subject finite difference.  Ground truth is the same model written as ODEs.
  .tvData <- function() {
    d <- nlmixr2data::theo_sd
    d$WT <- d$WT * exp(0.01 * d$TIME) # changes within each subject
    d
  }
  .tvModel <- function(lin) {
    if (lin) {
      function() {
        ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; wka <- 0.3; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
        model({
          ka <- exp(tka) * (WT / 70)^wka; cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
          cp <- linCmt(); cp ~ add(add.sd)
        })
      }
    } else {
      function() {
        ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; wka <- 0.3; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
        model({
          ka <- exp(tka) * (WT / 70)^wka; cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
          d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
          cp <- center / v; cp ~ add(add.sd)
        })
      }
    }
  }
  .tvCtl <- function(...) {
    foceiControl(
      print = 0L,
      covMethod = "",
      sigdig = 4,
      maxOuterIterations = 0L,
      maxInnerIterations = 300L,
      fast = TRUE,
      ...
    )
  }

  test_that("time-varying linCmt() parameter: carried outer gradient matches the ODE twin", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    skip_if_not(.foceiLinCmtHasSecondOrder())
    d <- .tvData()
    fl <- suppressMessages(suppressWarnings(nlmixr2(.tvModel(TRUE), d, "focei", .tvCtl())))
    fo <- suppressMessages(suppressWarnings(nlmixr2(
      .tvModel(FALSE),
      d,
      "focei",
      .tvCtl(rxControl = rxode2::rxControl(atol = 1e-10, rtol = 1e-10))
    )))
    expect_equal(fl$objf, fo$objf, tolerance = 1e-6)
    gl <- .foceiGradDirect(fl)
    go <- .foceiGradDirect(fo)
    expect_false(is.null(gl))
    expect_gt(fl$env$nAnalyticGradDirect, 0)
    expect_equal(fl$env$nOuterFdInd, 0L)
    nm <- names(fixef(fl))
    # tka/wka sit on the time-varying ka slot (kernel row 2 * ncmt)
    expect_true(all(abs(gl[c("tka", "wka")]) > 1e-3))
    expect_equal(unname(gl[nm]), unname(go[nm]), tolerance = 1e-6)
  })
})
