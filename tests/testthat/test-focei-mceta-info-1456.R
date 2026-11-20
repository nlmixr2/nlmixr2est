nmTest({
  .mcetaInfoMod <- function() {
    ini({
      tka <- 0.45
      tcl <- log(c(0, 2.7, 100))
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

  .mcetaInfoFit <- function(...) {
    suppressWarnings(suppressMessages(
      nlmixr2(
        .mcetaInfoMod,
        theo_sd,
        "focei",
        foceiControl(
          print = 0L,
          covMethod = "",
          calcTables = FALSE,
          maxOuterIterations = 2L,
          ...
        )
      )
    ))
  }

  .cols <- c(
    "startZero",
    "startDraw",
    "wonZero",
    "wonDraw",
    "ranked",
    "flipped",
    "time"
  )

  test_that("mceta <= 0 reports no multistart in fit$mcetaInfo (#1456)", {
    skip_on_cran()
    for (.m in c(0L, -2L)) {
      .f <- .mcetaInfoFit(mceta = .m)
      .info <- .f$mcetaInfo
      expect_equal(names(.info$total), .cols)
      expect_equal(.info$total[["startDraw"]], 0)
      expect_equal(.info$total[["wonDraw"]], 0)
      expect_equal(.info$total[["startZero"]], 0)
      expect_true(nrow(.info$iter) > 0)
      expect_true(all(.info$iter$startDraw == 0L))
    }
  })

  test_that("mceta = 5, mcetaOpt = TRUE counts every start (#1456)", {
    skip_on_cran()
    .nsub <- length(unique(theo_sd$ID))
    .f <- .mcetaInfoFit(mceta = 5L, mcetaOpt = TRUE)
    .info <- .f$mcetaInfo
    .it <- .info$iter
    expect_equal(.info$mceta, 5L)
    expect_true(.info$mcetaOpt)
    expect_equal(names(.it), c("iter", .cols))
    # a row holds whole searches (unrecorded probes fold into the next row)
    expect_true(all(.it$startZero %% .nsub == 0L))
    expect_true(all(.it$startDraw == 4L * .it$startZero))
    expect_true(all(.it$wonZero + .it$wonDraw == .it$startZero))
    expect_equal(.it$startZero[1], .nsub)
    for (.c in .cols) {
      expect_true(all(.it[[.c]] >= 0))
    }
    expect_true(all(.it$ranked >= .it$flipped))
    expect_true(.info$total[["ranked"]] >= .info$total[["flipped"]])
    expect_true(all(diff(.it$iter) > 0))
    # totals cover the history and agree with the kept-start counter
    expect_true(all(.info$total[.cols] >= colSums(.it[, .cols])))
    expect_equal(
      unname(.info$total[c("wonZero", "wonDraw")]),
      unname(as.numeric(.f$env$nMcetaStart[c("zero", "sample")]))
    )
  })

  test_that("mceta = 5 without mcetaOpt runs one start per subject (#1456)", {
    skip_on_cran()
    .nsub <- length(unique(theo_sd$ID))
    .it <- .mcetaInfoFit(mceta = 5L)$mcetaInfo$iter
    expect_true(all((.it$startZero + .it$startDraw) %% .nsub == 0L))
    expect_equal(.it$startZero, .it$wonZero)
    expect_equal(.it$startDraw, .it$wonDraw)
  })
})
