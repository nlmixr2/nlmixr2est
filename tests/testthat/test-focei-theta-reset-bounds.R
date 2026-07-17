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
    ctl <- foceiControl(resetThetaP = 0.4, resetThetaCheckPer = 1,
                        print = 0, maxOuterIterations = 40L,
                        covMethod = "", calcTables = FALSE)

    nReset <- 0L
    fit <- withCallingHandlers(
      suppressWarnings(nlmixr(boundedReset, theo_sd, est = "focei", control = ctl)),
      message = function(m) {
        if (grepl("ETA drift", conditionMessage(m), fixed = TRUE)) {
          nReset <<- nReset + 1L
        }
        invokeRestart("muffleMessage")
      })

    # the reset machinery must have run (otherwise this is not testing #454)
    expect_gt(nReset, 0L)

    idf <- fit$ui$iniDf
    th <- idf[!is.na(idf$ntheta), ]
    # every estimated population parameter must respect its bounds
    expect_true(all(th$est >= th$lower - 1e-6))
    expect_true(all(th$est <= th$upper + 1e-6))
    # the tightly-bounded parameter in particular stays at/under its upper bound
    expect_lte(th$est[th$name == "tcl"], 0.2 + 1e-6)
  })
})
