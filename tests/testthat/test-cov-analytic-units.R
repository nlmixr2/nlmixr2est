# Always-run units of the analytic covariance's observation rows and decline
# reasons (R/foceiCovAnalytic.R), none of which compiles a model, and the static
# gates on the decline chain (tests/testthat/helper-cov-analytic-declines.R); the
# fits are in test-cov-analytic-obs-rows.R, -reason.R and -declines.R (weekly).

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

  test_that("a decline with no reason gives the internal-error sentinel, and the fallback runs", {
    # unreachable while every decline records a reason; mocked to pin the wording
    local_mocked_bindings(.foceiCalcRanalyticAssemble = function(e) NULL)
    .e <- new.env(parent = emptyenv())
    expect_null(.foceiCalcRanalytic(.e))
    expect_identical(.e$.analyticReason, .foceiAnalyticNoReason)
    expect_identical(.e$.analyticReason, "an internal error (please report it)")
    expect_identical(.e$.analyticReasonShort, "an internal error (please report it)")
    expect_identical(
      .e$.analyticWarnNoFallback,
      "\"analytic\" not usable (an internal error (please report it)); no fallback"
    )
    expect_identical(.e$covTried, data.frame(method = "analytic", outcome = .foceiAnalyticNoReason))
    # setCov() stops with it
    local_mocked_bindings(.foceiCovAnalyticCalc = function(fit) NULL)
    .foceiAnalyticWhyReset()
    expect_error(
      .setCovAnalytic(NULL, new.env(parent = emptyenv()), "analytic"),
      "covMethod=\"analytic\" could not be computed for this fit (an internal error (please report it)); the covariance is left unchanged",
      fixed = TRUE
    )
  })

  test_that("every exit of the entry points and assemblers is a named decline", {
    .other <- function(fn, nested = FALSE) {
      .k <- unname(.covAnalyticReturnKinds(fn, nested = nested))
      sort(.k[.k != ".foceiAnalyticFallback"])
    }
    expect_identical(.other(".foceiCovAnalyticCalc"), character(0))
    # the covType opt-out, which the C++ hook never reaches
    expect_identical(.other(".foceiCalcRanalyticAssemble"), "NULL")
    # the floored-R signal, twice
    expect_identical(.other(".foceiAnalyticAssembleR"), c("structure", "structure"))
    # the FOCE branch's result
    expect_identical(.other(".foceiAnalyticAssembleRFR"), "R")
    # the one return() of a function they define: the quadrature-node callback's
    # NULL, which sets the floored flag its caller turns into that signal
    expect_identical(unname(.covAnalyticReturnKinds(".foceiAnalyticAssembleR", nested = TRUE)), "NULL")
    for (.fn in c(".foceiCovAnalyticCalc", ".foceiCalcRanalyticAssemble", ".foceiAnalyticAssembleRFR")) {
      expect_identical(.other(.fn, nested = TRUE), character(0), label = .fn)
    }
  })

  test_that("the closure walk reads calls and ui properties, not list fields", {
    .acc <- new.env(parent = emptyenv())
    .acc$syms <- character(0)
    .acc$fields <- character(0)
    .covAnalyticWalk(
      quote({
        .a <- ui$predDfFocei
        .b <- rxode2::rxUiGet(ui, "foce")
        f(x$g$h, y@k, z[, 1])
        do.call(".foceiX", list(1))
      }),
      .acc
    )
    expect_setequal(.acc$fields, c("predDfFocei", "foce", "g", "h"))
    # the fields after $ and @ are not names; a string is (do.call by name)
    expect_setequal(
      .acc$syms,
      c(
        "{",
        "<-",
        ".a",
        ".b",
        "::",
        "rxode2",
        "rxUiGet",
        "ui",
        "foce",
        "f",
        "x",
        "y",
        "[",
        "z",
        "do.call",
        ".foceiX",
        "list"
      )
    )
    .ns <- asNamespace("nlmixr2est")
    expect_true("rxUiGet.loadPruneSens" %in% .covAnalyticNamedFns(".foceiAnalyticAugModelDirs", .ns))
  })

  test_that("the closure walk stops at a leaf and follows a call by name", {
    .fake <- new.env(parent = emptyenv())
    .fake$entry <- function() {
      leafA()
      other()
    }
    .fake$leafA <- function() below()
    .fake$below <- function() NULL
    .fake$other <- function() do.call("byName", list())
    .fake$byName <- function() NULL
    expect_setequal(
      .covAnalyticDeclineClosure(.fake, entries = "entry", leaves = "leafA"),
      c("entry", "leafA", "other", "byName")
    )
    expect_setequal(
      .covAnalyticDeclineClosure(.fake, stopAtLeaves = FALSE, entries = "entry", leaves = "leafA"),
      c("entry", "leafA", "below", "other", "byName")
    )
  })

  test_that("every helper of the decline chain that can give up is a leaf or exempt for a reason", {
    .full <- .covAnalyticDeclineClosure(stopAtLeaves = FALSE)
    .reach <- .covAnalyticDeclineClosure()
    .leaves <- unique(unlist(.covAnalyticDeclineLeaves))
    expect_true(all(.leaves %in% .full), label = paste(setdiff(.leaves, .full), collapse = ", "))
    .gives <- .reach[vapply(.reach, .covAnalyticGivesUp, logical(1))]
    # a new helper that can give up must be made a leaf or exempted; a stale
    # exemption fails too
    expect_setequal(setdiff(.gives, .leaves), names(.covAnalyticDeclineExempt))
    # reached only through leaves: the solve's row check and the model build
    # behind .foceiAnalyticAugModelDirs
    .below <- setdiff(.full, .reach)
    .known <- c(
      ".foceiAnalyticObsAligned",
      "rxUiGet.loadPruneSens",
      "rxUiGet.foceiCmtPreModel",
      "rxUiGet.interpLinesStr"
    )
    expect_true(all(.known %in% .below), label = paste(setdiff(.known, .below), collapse = ", "))
    # fit$foceiControl is a list field, not a call of foceiControl()
    expect_false("foceiControl" %in% .full)
  })

  test_that("the pooled solve declines when its observation rows do not match", {
    .rows <- new.env(parent = emptyenv())
    .rows$n <- c(3L, 2L)
    local_mocked_bindings(
      .vaeOuterCols = function(am) list(),
      vaeOuterSolve_ = function(...) {
        structure(lapply(.rows$n, function(.n) list(a = matrix(0, .n, 1L))), ok = c(1L, 1L))
      }
    )
    .solve <- function() {
      .foceiAnalyticSolveAll(
        list(cores = 1L),
        c(THETA_1_ = 1),
        matrix(0, 2L, 1L),
        1:2,
        data.frame(),
        list(1:3, 1:3),
        1e-6,
        labels = c(101, 102)
      )
    }
    .foceiAnalyticWhyReset()
    expect_null(.solve())
    expect_identical(.foceiAnalyticWhy$detail, "subject 102: 2 observation rows solved for 3 observations")
    .rows$n <- c(3L, 3L)
    expect_length(.solve(), 2L)
  })

  # A per-subject solve whose A is linear in ETA_1_ (slope 0.5), failing where
  # `fails(eta)` is TRUE; `calls$x` records the eta of every solve
  .unitsSolveFA <- function(calls, fails) {
    function(aug, params, ev, times, tol = 1e-10) {
      .eta <- params[["ETA_1_"]]
      calls$x <- c(calls$x, .eta)
      if (fails(.eta)) {
        return(.foceiAnalyticNote("a test solve failure"))
      }
      list(f = c(1, 2), a = matrix(c(0.1, 0.2), 2L, 1L), A = array(0.01 + 0.5 * .eta, c(2L, 1L, 1L)))
    }
  }

  test_that("the third-order differencing declines when a side of its step fails", {
    # every solve but the base one fails, then only the plus side, then only the
    # minus side: shi21Central gives a zero, a zero and a forward difference
    .fails <- list(
      all = function(.e) .e != 0,
      plus = function(.e) .e > 0,
      minus = function(.e) .e < 0
    )
    for (.nm in names(.fails)) {
      .calls <- new.env(parent = emptyenv())
      .calls$x <- numeric(0)
      local_mocked_bindings(.foceiAnalyticSolveFA = .unitsSolveFA(.calls, .fails[[.nm]]))
      .foceiAnalyticWhyReset()
      expect_null(
        .foceiAnalyticSolveSubjectFD3(list(dirs = "ETA_1_"), c(ETA_1_ = 0), data.frame(), c(0, 1)),
        label = .nm
      )
      expect_identical(.foceiAnalyticWhy$detail, "a test solve failure", label = .nm)
      # the check of the step's sides reuses the differencing's own solves
      expect_identical(anyDuplicated(.calls$x), 0L, label = .nm)
    }
  })

  test_that("the third-order differencing declines on a step it did not solve on both sides", {
    # shi21Central always solves the step it returns; a step it did not is a change
    # in it, so the subject declines instead of solving there
    .calls <- new.env(parent = emptyenv())
    .calls$x <- numeric(0)
    local_mocked_bindings(
      .foceiAnalyticSolveFA = .unitsSolveFA(.calls, function(.e) FALSE),
      shi21CentralWrap = function(f, t, f0, idx, ef) list(h = 0.125, gr = rep(0, length(f0)))
    )
    .foceiAnalyticWhyReset()
    expect_null(.foceiAnalyticSolveSubjectFD3(list(dirs = "ETA_1_"), c(ETA_1_ = 0), data.frame(), c(0, 1)))
    expect_identical(.foceiAnalyticWhy$detail, "a third-order step with a side the differencing did not solve")
    expect_identical(.calls$x, 0)
  })

  test_that("a differencing solve of another size counts as a failed one", {
    # shiWrapCall_ turns such a result into NaN; the side check must agree
    local_mocked_bindings(.foceiAnalyticSolveFA = function(aug, params, ev, times, tol = 1e-10) {
      .eta <- params[["ETA_1_"]]
      .n <- if (.eta == 0) 2L else 3L
      list(f = seq_len(.n), a = matrix(0.1, .n, 1L), A = array(0.01 + 0.5 * .eta, c(.n, 1L, 1L)))
    })
    .foceiAnalyticWhyReset()
    expect_null(.foceiAnalyticSolveSubjectFD3(list(dirs = "ETA_1_"), c(ETA_1_ = 0), data.frame(), c(0, 1)))
    expect_identical(.foceiAnalyticWhy$detail, "a non-finite third-order sensitivity")
  })

  test_that("a successful third-order differencing solves no point twice", {
    .calls <- new.env(parent = emptyenv())
    .calls$x <- numeric(0)
    local_mocked_bindings(.foceiAnalyticSolveFA = .unitsSolveFA(.calls, function(.e) FALSE))
    .foceiAnalyticWhyReset()
    .r <- .foceiAnalyticSolveSubjectFD3(list(dirs = "ETA_1_"), c(ETA_1_ = 0), data.frame(), c(0, 1))
    # A is linear in the eta, so the central difference is exact up to rounding
    expect_equal(.r$Ath, array(0.5, c(2L, 1L, 1L, 1L)), tolerance = 1e-12)
    expect_identical(anyDuplicated(.calls$x), 0L)
    expect_null(.foceiAnalyticWhy$detail)
  })
})
