nmTest({
  # A time-varying covariate on a linCmt() parameter: the outer gradient's columns come
  # from rxode2lincmt's direction carry (.foceiOuterLinCmtCarry) instead of the
  # parameter-space reads, which read NA there and sent every subject to the
  # per-subject finite difference.  Ground truth is the same model written as ODEs.
  .tvData <- function(d = nlmixr2data::theo_sd) {
    d$WT <- d$WT * exp(0.01 * d$TIME) # changes within each subject
    d
  }
  .tvModel <- function(
    lin,
    kaLine = "ka <- exp(tka) * (WT / 70)^wka",
    err = "add(add.sd)",
    errIni = "add.sd <- 0.7"
  ) {
    f <- eval(str2lang(sprintf(
      paste(
        "function() {",
        "ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; wka <- 0.3; %s eta.cl ~ 0.3; eta.v ~ 0.1; %s })",
        "model({ %s; cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)",
        "d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center",
        "cp <- center / v; cp ~ %s }) }",
        sep = "\n"
      ),
      if (grepl("eta.ka", kaLine, fixed = TRUE)) "eta.ka ~ 0.2;" else "",
      errIni,
      kaLine,
      err
    )))
    if (!lin) {
      return(f)
    }
    suppressMessages(rxode2::model(f, -d/dt(depot), -d/dt(center), cp <- linCmt()))
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
  # the linCmt() fit and its ODE twin: same objective and the same gradient, omegas included
  .tvCompare <- function(d, ...) {
    fl <- suppressMessages(suppressWarnings(nlmixr2(.tvModel(TRUE, ...), d, "focei", .tvCtl())))
    fo <- suppressMessages(suppressWarnings(nlmixr2(
      .tvModel(FALSE, ...),
      d,
      "focei",
      .tvCtl(rxControl = rxode2::rxControl(atol = 1e-10, rtol = 1e-10))
    )))
    expect_equal(fl$objf, fo$objf, tolerance = 1e-6)
    gl <- .foceiGradDirect(fl)
    go <- .foceiGradDirect(fo)
    expect_false(is.null(gl))
    nm <- intersect(names(gl), names(go))
    expect_gt(length(nm), length(fixef(fl)))
    expect_equal(unname(gl[nm]), unname(go[nm]), tolerance = 1e-6)
    list(fl = fl, gl = gl)
  }

  test_that("time-varying linCmt() parameter: carried outer gradient matches the ODE twin", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    skip_if_not(.foceiLinCmtHasSecondOrder())
    r <- .tvCompare(.tvData())
    expect_gt(r$fl$env$nAnalyticGradDirect, 0)
    expect_equal(r$fl$env$nOuterFdInd, 0L)
    # tka/wka sit on the time-varying ka slot (kernel row 2 * ncmt)
    expect_true(all(abs(r$gl[c("tka", "wka")]) > 1e-3))
  })

  test_that("carried outer gradient with multiple doses and proportional error", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    skip_if_not(.foceiLinCmtHasSecondOrder())
    d <- .tvData(nlmixr2data::theo_md)
    r <- .tvCompare(d, err = "prop(prop.sd)", errIni = "prop.sd <- 0.2")
    expect_gt(r$fl$env$nAnalyticGradDirect, 0)
    expect_equal(r$fl$env$nOuterFdInd, 0L)
  })

  test_that("an eta on the time-varying slot is declined, and the gradient stays right", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    skip_if_not(.foceiLinCmtHasSecondOrder())
    .ka <- "ka <- exp(tka + eta.ka) * (WT / 70)^wka"
    .ui <- rxode2::rxUiDecompress(.tvModel(TRUE, kaLine = .ka)())
    rxode2::rxAssignControlValue(.ui, "foceiConstCovs", character(0))
    .s <- .ui$loadPruneSens
    expect_equal(
      .foceiOuterLinCmtCarry(.ui, .s, c("ETA_1_", "ETA_2_", "ETA_3_", "THETA_4_")),
      "an eta on a time-varying linCmt() parameter"
    )
    # the same slot without its eta is carried
    .ui <- rxode2::rxUiDecompress(.tvModel(TRUE)())
    rxode2::rxAssignControlValue(.ui, "foceiConstCovs", character(0))
    expect_true(is.list(.foceiOuterLinCmtCarry(.ui, .ui$loadPruneSens, c("ETA_1_", "ETA_2_", "THETA_4_"))))
    d <- .tvData()
    fl <- suppressMessages(suppressWarnings(nlmixr2(.tvModel(TRUE, kaLine = .ka), d, "focei", .tvCtl())))
    fo <- suppressMessages(suppressWarnings(nlmixr2(
      .tvModel(FALSE, kaLine = .ka),
      d,
      "focei",
      .tvCtl(rxControl = rxode2::rxControl(atol = 1e-10, rtol = 1e-10))
    )))
    expect_equal(fl$objf, fo$objf, tolerance = 1e-6)
  })
})
