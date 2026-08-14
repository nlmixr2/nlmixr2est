nmTest({
  # #923: an over-parameterized saem run can end with a singular (all-zero,
  # non-finite or negative-definite) Omega.  nmNearPD() itself errors on those,
  # which aborted the completed fit at the post-processing step.  The fit now
  # degrades (with a $runInfo note) instead of erroring.

  .om <- function(x, nm=c("eta.ka", "eta.cl")) {
    matrix(x, 2, 2, dimnames=list(nm, nm))
  }

  test_that(".foceiRepairOmega always returns a cholesky-able omega", {
    # the fully degenerate cases: nmNearPD() errors on all three
    .want <- c("singular omega", "non-finite", "singular omega")
    .degs <- list(.om(0), .om(c(0.4, NaN, NaN, 0.2)), .om(c(-1, 0, 0, -2)))
    for (.i in seq_along(.degs)) {
      .deg <- .degs[[.i]]
      expect_error(nmNearPD(.deg))
      expect_warning(.r <- .foceiRepairOmega(.deg), .want[.i])
      expect_false(inherits(try(chol(.r), silent=TRUE), "try-error"))
      expect_equal(dimnames(.r), dimnames(.deg))
    }
    # a merely singular omega is still repaired by nearPD (no fallback warning)
    expect_warning(.r <- .foceiRepairOmega(.om(c(0.5, 0, 0, 0))), NA)
    expect_false(inherits(try(chol(.r), silent=TRUE), "try-error"))
    expect_equal(dimnames(.r), dimnames(.om(0)))
    # a good omega is returned unchanged
    expect_equal(.foceiRepairOmega(.om(c(0.5, 0, 0, 0.2))), .om(c(0.5, 0, 0, 0.2)))
  })

  test_that("the focei post-processing env survives an all-zero omega (#923)", {
    .mod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.mod))
    # what a collapsed saem Omega looks like coming out of .getSaemOmega()
    .iniDf <- .ui$iniDf
    .iniDf$est[!is.na(.iniDf$neta1)] <- 0
    assign("iniDf", .iniDf, envir=.ui)
    assign("control", foceiControl(), envir=.ui)
    expect_error(nmNearPD(.ui$omega))

    .env <- new.env(parent=emptyenv())
    .env$etaNames <- c("eta.ka", "eta.cl")
    expect_warning(.foceiOptEnvSetupBounds(.ui, .env), "singular omega")
    # the mechanism that used to abort the fit now produced a usable inverse
    expect_false(is.null(.env$rxInv))
  })

  test_that(".saemWarnDegenerateOmega notes a collapsed omega in $runInfo", {
    expect_warning(.saemWarnDegenerateOmega(.om(0)), "all omega variances are zero")
    expect_warning(.saemWarnDegenerateOmega(.om(c(0.5, 0, 0, 0))),
                   "omega variance collapsed to zero: eta.cl")
    expect_warning(.saemWarnDegenerateOmega(.om(c(0.5, NaN, NaN, 0.2))), "non-finite")
    # a user-fixed zero variance is a modeling choice, not a collapse
    expect_warning(.saemWarnDegenerateOmega(.om(c(0.5, 0, 0, 0)), fixed="eta.cl"), NA)
    expect_warning(.saemWarnDegenerateOmega(.om(c(0.5, 0, 0, 0.2))), NA)
    expect_warning(.saemWarnDegenerateOmega(matrix(numeric(0), 0, 0)), NA)
  })

  test_that(".saemCreateOutput degrades to a table-free fit", {
    .calls <- new.env(parent=emptyenv())
    .calls$n <- 0L
    .calls$calcTables <- logical(0)
    .mockBuild <- function(nFail) {
      function(ui, data, control, table, env, est) {
        .calls$n <- .calls$n + 1L
        .calls$calcTables <- c(.calls$calcTables, control$calcTables)
        env$leftBehind <- TRUE
        if (.calls$n <= nFail) stop("build failure ", .calls$n)
        "fit"
      }
    }
    .mkEnv <- function() {
      .e <- new.env(parent=emptyenv())
      .e$ui <- NULL
      .e$origData <- NULL
      .e$control <- list(calcTables=TRUE)
      .e$table <- list(cwres=TRUE, npde=TRUE)
      .e
    }

    # first attempt fails -> retried without tables, warns, keeps the fit
    local_mocked_bindings(nlmixr2CreateOutputFromUi=.mockBuild(1L))
    .env <- .mkEnv()
    expect_warning(.ret <- .saemCreateOutput(.env), "fit degraded, no tables")
    expect_equal(.ret, "fit")
    expect_equal(.calls$n, 2L)
    # the retry really did turn the table step off ...
    expect_equal(.calls$calcTables, c(TRUE, FALSE))
    expect_false(.env$table$cwres)
    expect_false(.env$table$npde)

    # ... and the failed attempt's writes into the environment were rolled back
    # before the retry (`leftBehind` only comes from the successful second call)
    .calls$n <- 0L
    .calls$calcTables <- logical(0)
    .env2 <- .mkEnv()
    local_mocked_bindings(nlmixr2CreateOutputFromUi=function(ui, data, control, table, env, est) {
      .calls$n <- .calls$n + 1L
      if (.calls$n == 1L) {
        env$leftBehind <- TRUE
        stop("build failure")
      }
      exists("leftBehind", envir=env, inherits=FALSE)
    })
    expect_warning(.ret2 <- .saemCreateOutput(.env2), "fit degraded")
    expect_false(.ret2)

    # both attempts fail -> the original error is what the user sees
    .calls$n <- 0L
    local_mocked_bindings(nlmixr2CreateOutputFromUi=.mockBuild(2L))
    expect_error(.saemCreateOutput(.mkEnv()), "build failure 1")
  })
})
