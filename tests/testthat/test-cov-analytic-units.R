# Always-run units of the analytic covariance's observation rows and decline
# reasons (R/foceiCovAnalytic.R); the fits that use them are in
# test-cov-analytic-obs-rows.R and test-cov-analytic-reason.R (weekly batch 2).

nmTest({
  test_that(".foceiAnalyticNote() keeps a detail on one line, and NULL clears it", {
    .foceiAnalyticWhyReset()
    expect_null(.foceiAnalyticNote("first line\n  second\tline  "))
    expect_identical(.foceiAnalyticWhy$detail, "first line second line")
    .foceiAnalyticNote(NULL)
    expect_null(.foceiAnalyticWhy$detail)
  })

  test_that(".foceiAnalyticFallback() records the step and its detail", {
    .foceiAnalyticWhyReset()
    expect_message(
      expect_null(.foceiAnalyticFallback("a feature")),
      "covType=\"analytic\": a feature is out of analytic-covariance scope; using the finite-difference covariance instead",
      fixed = TRUE
    )
    expect_identical(.foceiAnalyticWhy$reason, "a feature")
    expect_identical(.foceiAnalyticWhy$short, "a feature")

    # a failure appends what the failing solve noted, and consumes it
    .foceiAnalyticWhyReset()
    .foceiAnalyticNote("what the solve said")
    expect_message(
      .foceiAnalyticFallback("a failed step", scope = FALSE),
      "covType=\"analytic\": a failed step (what the solve said); using the finite-difference covariance instead",
      fixed = TRUE
    )
    expect_identical(.foceiAnalyticWhy$reason, "a failed step (what the solve said)")
    expect_identical(.foceiAnalyticWhy$short, "a failed step")
    expect_null(.foceiAnalyticWhy$detail)

    # a scope reason does not take a detail
    .foceiAnalyticWhyReset()
    .foceiAnalyticNote("unrelated")
    suppressMessages(.foceiAnalyticFallback("a feature"))
    expect_identical(.foceiAnalyticWhy$reason, "a feature")

    # keep = TRUE keeps a reason a callee gave, without a second message
    .foceiAnalyticWhyReset()
    suppressMessages(.foceiAnalyticFallback("the callee's reason"))
    expect_silent(.foceiAnalyticFallback("the caller's reason", keep = TRUE))
    expect_identical(.foceiAnalyticWhy$reason, "the callee's reason")
    .foceiAnalyticWhyReset()
    expect_message(.foceiAnalyticFallback("the caller's reason", keep = TRUE), "the caller's reason", fixed = TRUE)
  })

  test_that(".foceiAnalyticClip() shortens to the width, ellipsis included", {
    expect_identical(.foceiAnalyticClip("short"), "short")
    .x <- strrep("x", 60)
    expect_identical(.foceiAnalyticClip(.x), paste0(strrep("x", 51), "..."))
    expect_identical(nchar(.foceiAnalyticClip(.x)), 54L)
    expect_identical(.foceiAnalyticClip(.x, 10L), "xxxxxxx...")
  })

  test_that(".foceiAnalyticObsAligned() accepts only the expected rows, and says why", {
    .foceiAnalyticWhyReset()
    expect_true(.foceiAnalyticObsAligned(c(0, 1, 2), c(0, 1, 2)))
    expect_null(.foceiAnalyticWhy$detail)
    expect_false(.foceiAnalyticObsAligned(numeric(0), c(0, 1)))
    expect_identical(.foceiAnalyticWhy$detail, "the augmented solve returned no observation rows")
    expect_false(.foceiAnalyticObsAligned(c(0, 1, 1, 2), c(0, 1, 2)))
    expect_identical(.foceiAnalyticWhy$detail, "4 observation rows solved for 3 observations")
    expect_false(.foceiAnalyticObsAligned(c(0, 1, 2.5), c(0, 1, 2)))
    expect_identical(.foceiAnalyticWhy$detail, "solved observation times differ from the data's")
  })

  test_that("the subject solve keeps the observation records, whatever else shares their time", {
    skip_on_cran()
    .mod <- function() {
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
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::assertRxUi(.mod))
    .am <- .foceiAnalyticAugModelDirs(.ui, paste0("ETA_", 1:3, "_"))
    expect_false(is.null(.am))
    .thr <- .ui$iniDf[!is.na(.ui$iniDf$ntheta), ]
    .thr <- .thr[order(.thr$ntheta), ]
    .p <- c(
      stats::setNames(.thr$est, paste0("THETA_", .thr$ntheta, "_")),
      stats::setNames(rep(0, 3), paste0("ETA_", 1:3, "_"))
    )
    .s1 <- nlmixr2data::theo_sd[nlmixr2data::theo_sd$ID == 1, ]
    .s1$MDV <- as.integer(.s1$EVID != 0)
    .t <- .s1$TIME[.s1$EVID == 0]
    .add <- function(evid, dv, mdv) {
      .r <- .s1[.s1$EVID == 0, ][1, ]
      .r$TIME <- .t[3]
      .r$EVID <- evid
      .r$DV <- dv
      .r$MDV <- mdv
      .d <- rbind(.s1, .r)
      .d[order(.d$TIME, -.d$EVID), ]
    }
    .e0 <- .foceiAnalyticSolveFA(.am, .p, .s1, .t)
    expect_length(.e0$f, length(.t))
    # the predictions agree to 1e-9: measured 0 for the missing-DV and MDV = 1
    # records and 1.1e-11 for the EVID 2 one, whose compartment event restarts the
    # solver at that time
    for (.d in list(.add(2L, NA_real_, 1L), .add(0L, NA_real_, 0L), .add(0L, 5, 1L))) {
      .e <- .foceiAnalyticSolveFA(.am, .p, .d, .t)
      expect_length(.e$f, length(.t))
      expect_equal(.e$f, .e0$f, tolerance = 1e-9)
    }
    # the times only confirm the rows; a mismatch is a NULL with the mismatch noted
    .foceiAnalyticWhyReset()
    expect_null(.foceiAnalyticSolveFA(.am, .p, .s1, c(.t, 30)))
    expect_identical(
      .foceiAnalyticWhy$detail,
      sprintf("%d observation rows solved for %d observations", length(.t), length(.t) + 1L)
    )
    # a failed rxSolve() is noted on one line
    .foceiAnalyticWhyReset()
    expect_null(.foceiAnalyticRxSolve(.am, .p[-1], .s1, 1e-8))
    expect_match(.foceiAnalyticWhy$detail, "^rxSolve: ")
    expect_false(grepl("\n", .foceiAnalyticWhy$detail, fixed = TRUE))
  })
})
