nmTest({
  # Regression test for the "for subject(s): Unknown" bug.  rxode2 aggregates
  # ODE-solve warnings (lsoda nhnil / dop853 / intdy window-miss) and labels
  # them by subject via rxGetId(), which reads rxode2's global ID factor table.
  # FOCEi/SAEM pass a declassed data.frame to rxSolve_, so rxode2 never
  # populated that table during estimation and every id resolved to the literal
  # "Unknown".  The fix pushes the fit's idLvl into rxode2 (rxSetIdLvlFactors)
  # right after estimation setup (src/inner.cpp, src/saem.cpp).
  #
  # These tests use two internal rxode2 test entry points (added alongside the
  # fix) to assert the factor table is populated with the real ids after a
  # fit: _rxTestSolveWarnLabels clears/sets the table, _rxTestGetIdLabels reads
  # rxGetId() back.  They are skipped when run against an rxode2 that predates
  # those helpers (the matching rxode2 ships them).

  .haveProbe <- tryCatch({
    invisible(.Call("_rxTestGetIdLabels", 0L, PACKAGE = "rxode2"))
    TRUE
  }, error = function(e) FALSE)

  .clearFactors <- function() {
    invisible(.Call("_rxTestSolveWarnLabels", character(0), integer(0),
                    PACKAGE = "rxode2"))
  }
  .idLabels <- function(ids) {
    .Call("_rxTestGetIdLabels", as.integer(ids), PACKAGE = "rxode2")
  }

  .mkModel <- function() {
    function() {
      ini({
        lka <- 0.5
        lcl <- -3.2
        lv <- -1
        prop.err <- 0.1
        eta.cl ~ 0.1
      })
      model({
        ka <- exp(lka)
        cl <- exp(lcl + eta.cl)
        v <- exp(lv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - (cl / v) * center
        cp <- center / v
        cp ~ prop(prop.err)
      })
    }
  }
  # Known, non-sequential ids so a real label can't be confused with a 0-based
  # internal index.
  .mkData <- function() {
    d <- nlmixr2data::theo_sd
    d <- d[d$ID %in% 1:3, ]
    d$ID <- as.character(factor(d$ID, labels = c("101", "202", "303")))
    d
  }

  test_that("focei estimation populates rxode2's subject-id factors (not 'Unknown')", {
    skip_if_not(.haveProbe, "rxode2 lacks the _rxTest* warn helpers")
    d <- .mkData()
    # Clear the global factor table so only estimation can repopulate it; with
    # calcTables = FALSE there is no post-fit solve to mask the estimation
    # state.  Without the setup-time rxSetIdLvlFactors call these stay "Unknown".
    .clearFactors()
    expect_equal(.idLabels(0L), "Unknown")
    suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(
        .mkModel(), d, est = "focei",
        control = nlmixr2est::foceiControl(
          print = 0, maxOuterIterations = 2, maxInnerIterations = 5,
          calcTables = FALSE, covMethod = ""
        )
      )
    ))
    expect_equal(.idLabels(c(0L, 1L, 2L)), c("101", "202", "303"))
  })

  test_that("saem estimation populates rxode2's subject-id factors (not 'Unknown')", {
    skip_if_not(.haveProbe, "rxode2 lacks the _rxTest* warn helpers")
    d <- .mkData()
    .clearFactors()
    expect_equal(.idLabels(0L), "Unknown")
    suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(
        .mkModel(), d, est = "saem",
        control = nlmixr2est::saemControl(
          print = 0, nBurn = 10, nEm = 10, calcTables = FALSE
        )
      )
    ))
    expect_equal(.idLabels(c(0L, 1L, 2L)), c("101", "202", "303"))
  })
})
