nmTest({
  .scModel <- function() {
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
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  test_that("an analytic precursor seeds the full stage and the shortcut skips its off-diagonals", {
    skip_on_cran()
    # one analytic fit per request, so no request sees another's stored result
    .ctl <- foceiControl(print = 0, calcTables = FALSE, covMethod = "analytic")
    .none <- suppressMessages(.nlmixr(.scModel, theo_sd, "focei", .ctl))
    .sc <- suppressMessages(.nlmixr(.scModel, theo_sd, "focei", .ctl))
    .off <- suppressMessages(.nlmixr(.scModel, theo_sd, "focei", .ctl))
    suppressMessages(suppressWarnings(setCov(.none, "r,s (full)", rsControl(covPrecursor = NULL))))
    suppressMessages(suppressWarnings(setCov(.sc, "r,s (full)", rsControl(covPrecursor = "analytic"))))
    suppressMessages(suppressWarnings(
      setCov(.off, "r,s (full)", rsControl(covPrecursor = "analytic", covShortcut = FALSE))
    ))
    for (.f in list(.none, .sc, .off)) {
      expect_identical(.f$covMethod, "r,s (full)")
    }
    expect_null(.none$env$covPrecursorUsed[["r,s (full)"]])
    .rec <- .sc$env$covPrecursorUsed[["r,s (full)"]]
    expect_identical(.rec[c("source", "shortcut")], list(source = "analytic", shortcut = "accepted"))
    # four check directions, each explained to 1%
    expect_identical(dim(.rec$checks), c(4L, 3L))
    expect_true(all(abs(.rec$checks[, "measured"] - .rec$checks[, "predicted"]) <= .rec$checks[, "allowance"]))
    expect_equal(unname(.rec$checks[, "allowance"]), 0.01 * abs(unname(.rec$checks[, "predicted"])))
    expect_identical(
      .off$env$covPrecursorUsed[["r,s (full)"]][c("source", "shortcut")],
      list(source = "analytic", shortcut = "off")
    )
    # the predicted off-diagonals give the measured ones' covariance to within the
    # finite-difference error (0.4% on theo_sd)
    expect_equal(sqrt(diag(.sc$cov)), sqrt(diag(.none$cov)), tolerance = 0.01)
    expect_equal(sqrt(diag(.off$cov)), sqrt(diag(.none$cov)), tolerance = 0.01)
    expect_output(print(.sc), "from the \"analytic\" precursor (seeded steps; shortcut accepted)", fixed = TRUE)
  })

  test_that("the shortcut falls back to measuring the off-diagonals a precursor does not explain", {
    skip_on_cran()
    # the analytic covariance with its correlations removed predicts no off-diagonal
    # curvature, which the checks see; falling back measures them at the same steps,
    # so the covariance is the one computed without the shortcut
    .ctl <- foceiControl(print = 0, calcTables = FALSE, covMethod = "analytic")
    .sc <- suppressMessages(.nlmixr(.scModel, theo_sd, "focei", .ctl))
    .off <- suppressMessages(.nlmixr(.scModel, theo_sd, "focei", .ctl))
    for (.f in list(.sc, .off)) {
      .d <- diag(diag(.f$cov))
      dimnames(.d) <- dimnames(.f$cov)
      assign("cov", .d, envir = .f$env)
    }
    suppressMessages(suppressWarnings(setCov(.sc, "r,s (full)", rsControl(covPrecursor = "analytic"))))
    suppressMessages(suppressWarnings(
      setCov(.off, "r,s (full)", rsControl(covPrecursor = "analytic", covShortcut = FALSE))
    ))
    .rec <- .sc$env$covPrecursorUsed[["r,s (full)"]]
    expect_identical(.rec[c("source", "shortcut")], list(source = "analytic", shortcut = "fell back"))
    .last <- .rec$checks[nrow(.rec$checks), ]
    expect_gt(abs(.last[["measured"]] - .last[["predicted"]]), .last[["allowance"]])
    expect_identical(.sc$covMethod, "r,s (full)")
    expect_equal(.sc$cov, .off$cov, tolerance = 1e-10)
  })
})
