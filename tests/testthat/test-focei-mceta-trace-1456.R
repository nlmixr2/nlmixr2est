nmTest({
  .mcetaTraceMod <- function() {
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

  .mcetaTraceRun <- function(...) {
    .fit <- NULL
    .out <- utils::capture.output(
      # no suppressMessages(): it also silences the C-level iteration print
      .fit <- suppressWarnings(
        nlmixr2(
          .mcetaTraceMod,
          theo_sd,
          "focei",
          foceiControl(
            print = 1L,
            covMethod = "",
            calcTables = FALSE,
            maxOuterIterations = 2L,
            ...
          )
        )
      )
    )
    list(fit = .fit, trace = grep("^mceta:", .out, value = TRUE))
  }

  test_that("mcetaTrace is validated", {
    expect_error(foceiControl(mcetaTrace = NA))
    expect_error(foceiControl(mcetaTrace = c(TRUE, FALSE)))
    expect_false(foceiControl()$mcetaTrace)
  })

  test_that("mcetaTrace prints nothing without a multistart (#1456)", {
    skip_on_cran()
    for (.m in c(0L, -2L)) {
      expect_length(.mcetaTraceRun(mceta = .m, mcetaTrace = TRUE)$trace, 0L)
    }
  })

  test_that("mcetaTrace prints counts and leaves the fit unchanged (#1456)", {
    skip_on_cran()
    .nsub <- length(unique(theo_sd$ID))
    .on <- .mcetaTraceRun(mceta = 5L, mcetaOpt = TRUE, mcetaTrace = TRUE)
    .off <- .mcetaTraceRun(mceta = 5L, mcetaOpt = TRUE)
    expect_length(.off$trace, 0L)
    expect_true(length(.on$trace) > 0L)
    .first <- regmatches(
      .on$trace[1],
      regexec("starts ([0-9]+) zero \\+ ([0-9]+) draw", .on$trace[1])
    )[[1]]
    expect_equal(as.integer(.first[3]), 4L * as.integer(.first[2]))
    expect_equal(as.integer(.first[2]) %% .nsub, 0L)
    expect_equal(.on$fit$objf, .off$fit$objf)
    expect_equal(.on$fit$theta, .off$fit$theta)
  })
})
