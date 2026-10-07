nmTest({
  # #923: an over-parameterized saem run can end with a singular (all-zero,
  # non-finite or negative-definite) Omega.  nmNearPD() itself errors on those,
  # which aborted the completed fit at the post-processing step.  The fit now
  # degrades (with a $runInfo note) instead of erroring.

  .om <- function(x, nm = c("eta.ka", "eta.cl")) {
    matrix(x, 2, 2, dimnames = list(nm, nm))
  }

  test_that(".foceiRepairOmega always returns a cholesky-able omega", {
    # the fully degenerate cases: nmNearPD() errors on all three
    .want <- c("singular omega", "non-finite", "singular omega")
    .degs <- list(.om(0), .om(c(0.4, NaN, NaN, 0.2)), .om(c(-1, 0, 0, -2)))
    for (.i in seq_along(.degs)) {
      .deg <- .degs[[.i]]
      expect_error(nmNearPD(.deg))
      expect_warning(.r <- .foceiRepairOmega(.deg), .want[.i])
      expect_false(inherits(try(chol(.r), silent = TRUE), "try-error"))
      expect_equal(dimnames(.r), dimnames(.deg))
    }
    # a merely singular omega is still repaired by nearPD, and says so
    expect_warning(.r <- .foceiRepairOmega(.om(c(0.5, 0, 0, 0))), "not positive definite")
    expect_false(inherits(try(chol(.r), silent = TRUE), "try-error"))
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
    assign("iniDf", .iniDf, envir = .ui)
    assign("control", foceiControl(), envir = .ui)
    expect_error(nmNearPD(.ui$omega))

    .env <- new.env(parent = emptyenv())
    .env$etaNames <- c("eta.ka", "eta.cl")
    expect_warning(.foceiOptEnvSetupBounds(.ui, .env), "singular omega")
    # the mechanism that used to abort the fit now produced a usable inverse
    expect_false(is.null(.env$rxInv))
  })

  test_that("a degenerate saem omega is reported as estimated, not as its repair (issue 1140)", {
    .mod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl + eta.v ~ c(0.3, 0.01, 0.1)
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    # saem's estimate, made singular: eta.cl and eta.v perfectly correlated
    .cap <- new.env(parent = emptyenv())
    .getSaemOmega0 <- .getSaemOmega
    local_mocked_bindings(.getSaemOmega = function(env) {
      .getSaemOmega0(env)
      .om <- env$omega
      .om["eta.cl", "eta.v"] <- .om["eta.v", "eta.cl"] <- sqrt(.om["eta.cl", "eta.cl"] * .om["eta.v", "eta.v"])
      env$omega <- .om
      .cap$omega <- .om
      invisible()
    })
    fit <- suppressMessages(suppressWarnings(nlmixr2(
      .mod,
      theo_sd,
      "saem",
      saemControl(print = 0, nBurn = 5, nEm = 5, seed = 42, covMethod = "")
    )))
    expect_true(inherits(try(chol(.cap$omega), silent = TRUE), "try-error"))
    expect_equal(fit$omega, .cap$omega)
    expect_equal(fit$ui$omega, .cap$omega)
    # the tables needed the repair, and the fit says so
    expect_true(any(grepl("not positive definite", fit$runInfo, fixed = TRUE)))
  })

  test_that("a posthoc fit corrects a nearly singular ini() omega, but not a far one", {
    .mk <- function(cv) {
      .f <- function() {
        ini({
          tka <- 0.45
          tcl <- 1
          tv <- 3.45
          eta.ka ~ 0.6
          add.sd <- 0.7
        })
        model({
          ka <- exp(tka + eta.ka)
          cl <- exp(tcl + eta.cl)
          v <- exp(tv + eta.v)
          linCmt() ~ add(add.sd)
        })
      }
      suppressMessages(eval(bquote(rxode2::ini(.f, eta.cl + eta.v ~ c(0.1, .(cv), 0.1)))))
    }
    .ctl <- foceiControl(print = 0L, maxOuterIterations = 0L, covMethod = "", calcTables = FALSE)
    # a correlation of exactly 1, as a rounded import gives: corrected, and the
    # correction reported
    fit <- suppressMessages(.nlmixr(.mk(0.1), theo_sd, "focei", .ctl))
    .om <- fit$omega[c("eta.cl", "eta.v"), c("eta.cl", "eta.v")]
    expect_false(inherits(try(chol(.om), silent = TRUE), "try-error"))
    expect_equal(.om[1, 2], 0.1)
    expect_equal(diag(.om), c(eta.cl = 0.1, eta.v = 0.1), tolerance = 1e-3)
    expect_equal(fit$ui$omega, fit$omega)
    expect_true(any(grepl("nearly singular", fit$runInfo, fixed = TRUE)))
    # a correlation of 2 is no rounding: reported as given
    fit <- suppressMessages(.nlmixr(.mk(0.2), theo_sd, "focei", .ctl))
    expect_equal(fit$omega[["eta.cl", "eta.v"]], 0.2)
    expect_equal(fit$omega[["eta.cl", "eta.cl"]], 0.1)
    expect_true(any(grepl("not positive definite", fit$runInfo, fixed = TRUE)))
  })

  test_that(".saemWarnDegenerateOmega notes a collapsed omega in $runInfo", {
    expect_warning(.saemWarnDegenerateOmega(.om(0)), "all omega variances are zero")
    expect_warning(.saemWarnDegenerateOmega(.om(c(0.5, 0, 0, 0))), "omega variance collapsed to zero: eta.cl")
    expect_warning(.saemWarnDegenerateOmega(.om(c(0.5, NaN, NaN, 0.2))), "non-finite")
    # a user-fixed zero variance is a modeling choice, not a collapse
    expect_warning(.saemWarnDegenerateOmega(.om(c(0.5, 0, 0, 0)), fixed = "eta.cl"), NA)
    expect_warning(.saemWarnDegenerateOmega(.om(c(0.5, 0, 0, 0.2))), NA)
    expect_warning(.saemWarnDegenerateOmega(matrix(numeric(0), 0, 0)), NA)
  })

  test_that(".saemCreateOutput degrades to a table-free fit", {
    .calls <- new.env(parent = emptyenv())
    .calls$n <- 0L
    .calls$calcTables <- logical(0)
    .mockBuild <- function(nFail) {
      function(ui, data, control, table, env, est) {
        .calls$n <- .calls$n + 1L
        .calls$calcTables <- c(.calls$calcTables, control$calcTables)
        env$leftBehind <- TRUE
        if (.calls$n <= nFail) {
          stop("build failure ", .calls$n)
        }
        "fit"
      }
    }
    .mkEnv <- function() {
      .e <- new.env(parent = emptyenv())
      .e$ui <- NULL
      .e$origData <- NULL
      .e$control <- list(calcTables = TRUE)
      .e$table <- list(cwres = TRUE, npde = TRUE)
      .e
    }

    # first attempt fails -> retried without tables, warns, keeps the fit
    local_mocked_bindings(nlmixr2CreateOutputFromUi = .mockBuild(1L))
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
    local_mocked_bindings(nlmixr2CreateOutputFromUi = function(ui, data, control, table, env, est) {
      .calls$n <- .calls$n + 1L
      if (.calls$n == 1L) {
        env$leftBehind <- TRUE
        stop("build failure")
      }
      exists("leftBehind", envir = env, inherits = FALSE)
    })
    expect_warning(.ret2 <- .saemCreateOutput(.env2), "fit degraded")
    expect_false(.ret2)

    # both attempts fail -> the original error is what the user sees
    .calls$n <- 0L
    local_mocked_bindings(nlmixr2CreateOutputFromUi = .mockBuild(2L))
    expect_error(.saemCreateOutput(.mkEnv()), "build failure 1")
  })
})
