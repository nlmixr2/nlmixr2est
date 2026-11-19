nmTest({
  # Issue #454: a theta reset (soft mu-reference shift + restart) must never
  # move a population parameter outside its declared bounds.  The reset now
  # clamps the shifted theta into its (margin-adjusted) bounds instead of
  # skipping the shift, so a reset-heavy fit of a tightly-bounded mu-referenced
  # parameter is guaranteed to finish in range.
  test_that("theta resets keep population parameters within their bounds (#454)", {
    boundedReset <- function() {
      ini({
        tka <- 0.45
        tcl <- c(-2, -1, 0.2)   # tight upper bound (0.2) on log-CL
        tv <- 3.45
        add.sd <- 0.7
        eta.ka ~ 0.6
        eta.cl ~ 0.5
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }

    # aggressive reset settings so the theta-reset path is actually exercised
    ctl <- foceiControl(
      resetThetaP = 0.4,
      resetThetaCheckPer = 1,
      print = 0,
      maxOuterIterations = 40L,
      covMethod = "",
      calcTables = FALSE
    )

    nReset <- 0L
    fit <- withCallingHandlers(
      suppressWarnings(nlmixr(boundedReset, theo_sd, est = "focei", control = ctl)),
      message = function(m) {
        if (grepl("ETA drift", conditionMessage(m), fixed = TRUE)) {
          nReset <<- nReset + 1L
        }
        invokeRestart("muffleMessage")
      }
    )

    # the reset machinery must have run (otherwise this is not testing #454)
    expect_gt(nReset, 0L)
    # tcl's optimum (about 1) is past its bound, so the drift returns once the
    # optimizer moves it inward; the reset that would put it back is skipped
    expect_true(any(grepl("reset of 'tcl' skipped: back at its bound", fit$runInfo, fixed = TRUE)))

    idf <- fit$ui$iniDf
    th <- idf[!is.na(idf$ntheta), ]
    # every estimated population parameter must respect its bounds
    expect_true(all(th$est >= th$lower - 1e-6))
    expect_true(all(th$est <= th$upper + 1e-6))
    # the tightly-bounded parameter in particular stays at/under its upper bound
    expect_lte(th$est[th$name == "tcl"], 0.2 + 1e-6)
  })

  test_that("a theta pinned at its bound does not fire the theta reset again", {
    skip_on_cran()
    # The first reset moves tcl to its upper bound.  The restart's first
    # evaluation finds eta.cl drifting again, but tcl cannot follow, so no
    # reset fires there, even though eta.ka could take a tiny shift.  The outer
    # optimizer evaluates only its starting point, so every reset after the
    # first would come from that evaluation.
    .pinnedReset <- function() {
      ini({
        tka <- 0.45
        tcl <- c(-2, -1, 0.2)
        tv <- 3.45
        add.sd <- 0.7
        eta.ka ~ 0.6
        eta.cl ~ 0.5
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .opt <- function(par, fn, gr, lower, upper, control, ...) {
      .v <- fn(par)
      list(x = par, value = .v, convergence = 0L, message = "")
    }
    .ctl <- foceiControl(
      resetThetaP = 0.4,
      resetThetaCheckPer = 1,
      print = 0,
      covMethod = "",
      calcTables = FALSE,
      outerOpt = .opt
    )
    .acc <- new.env(parent = emptyenv())
    .acc$msg <- character(0)
    .fit <- withCallingHandlers(
      suppressWarnings(nlmixr(.pinnedReset, theo_sd, est = "focei", control = .ctl)),
      message = function(m) {
        .acc$msg <- c(.acc$msg, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
    expect_equal(sum(grepl("ETA drift", .acc$msg, fixed = TRUE)), 1L)
    expect_equal(fixef(.fit)[["tcl"]], 0.2, tolerance = 1e-5)
  })
})
