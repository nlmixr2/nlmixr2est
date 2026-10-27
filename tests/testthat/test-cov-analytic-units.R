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
    # an empty string and one too long for exists() are not names
    .fake$strings <- function() NULL
    body(.fake$strings) <- call("c", "", strrep("x", 10001L))
    expect_identical(.covAnalyticNamedFns("strings", .fake), character(0))
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

  # the pooled population solve of two subjects with three observations each
  .unitsSolveAllTwo <- function() {
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

  test_that("the pooled solve declines when its observation rows do not match", {
    .rows <- new.env(parent = emptyenv())
    .rows$n <- c(3L, 2L)
    local_mocked_bindings(
      .vaeOuterCols = function(am) list(),
      vaeOuterSolve_ = function(...) {
        structure(lapply(.rows$n, function(.n) list(a = matrix(0, .n, 1L))), ok = c(1L, 1L))
      }
    )
    .foceiAnalyticWhyReset()
    expect_null(.unitsSolveAllTwo())
    expect_identical(.foceiAnalyticWhy$detail, "subject 102: 2 observation rows solved for 3 observations")
    .rows$n <- c(3L, 3L)
    expect_length(.unitsSolveAllTwo(), 2L)
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
    expect_identical(.foceiAnalyticWhy$detail, "a third-order solve of another size")
  })

  test_that("a right-sized differencing solve with a non-finite value counts as failed", {
    local_mocked_bindings(.foceiAnalyticSolveFA = function(aug, params, ev, times, tol = 1e-10) {
      .eta <- params[["ETA_1_"]]
      list(f = c(1, 2), a = matrix(c(0.1, 0.2), 2L, 1L), A = array(if (.eta == 0) 0.01 else NaN, c(2L, 1L, 1L)))
    })
    .foceiAnalyticWhyReset()
    expect_null(.foceiAnalyticSolveSubjectFD3(list(dirs = "ETA_1_"), c(ETA_1_ = 0), data.frame(), c(0, 1)))
    expect_identical(.foceiAnalyticWhy$detail, "a non-finite third-order sensitivity")
  })

  test_that(".foceiAnalyticFD3Flat() flattens A, and AR with withR", {
    .E <- list(A = array(1:4, c(1L, 2L, 2L)), AR = array(5:8, c(1L, 2L, 2L)))
    expect_identical(.foceiAnalyticFD3Flat(.E, FALSE), 1:4)
    expect_identical(.foceiAnalyticFD3Flat(.E, TRUE), 1:8)
    expect_identical(.foceiAnalyticFD3Flat(list(A = .E$A), TRUE), 1:4)
  })

  # A per-subject solve whose A and AR are linear in ETA_1_ (slopes 0.5 and 0.25);
  # with `dropAR$on`, a perturbed solve leaves AR out
  .unitsSolveAR <- function(dropAR) {
    function(aug, params, ev, times, tol = 1e-10) {
      .eta <- params[["ETA_1_"]]
      .r <- list(f = c(1, 2), a = matrix(c(0.1, 0.2), 2L, 1L), A = array(0.01 + 0.5 * .eta, c(2L, 1L, 1L)))
      if (.eta == 0 || !isTRUE(dropAR$on)) {
        .r$AR <- array(0.02 + 0.25 * .eta, c(2L, 1L, 1L))
      }
      .r
    }
  }

  test_that("the third-order differencing differences AR with withR", {
    # linear A and AR, so the differences are exact
    .dropAR <- new.env(parent = emptyenv())
    .dropAR$on <- FALSE
    local_mocked_bindings(.foceiAnalyticSolveFA = .unitsSolveAR(.dropAR))
    .foceiAnalyticWhyReset()
    .r <- .foceiAnalyticSolveSubjectFD3(list(dirs = "ETA_1_"), c(ETA_1_ = 0), data.frame(), c(0, 1), withR = TRUE)
    expect_equal(.r$Ath, array(0.5, c(2L, 1L, 1L, 1L)), tolerance = 1e-12)
    expect_equal(.r$AthR, array(0.25, c(2L, 1L, 1L, 1L)), tolerance = 1e-12)
    # a perturbed solve without AR is shorter than f0, so it is one of another size
    .dropAR$on <- TRUE
    .foceiAnalyticWhyReset()
    expect_null(.foceiAnalyticSolveSubjectFD3(
      list(dirs = "ETA_1_"),
      c(ETA_1_ = 0),
      data.frame(),
      c(0, 1),
      withR = TRUE
    ))
    expect_identical(.foceiAnalyticWhy$detail, "a third-order solve of another size")
  })

  # A per-subject solve with one residual sigma (rx_rsig_): every slot affine in
  # ETA_1_, so the differences are exact.  Two observations; the values at eta 0
  # are what the expanded slots must hold.
  .unitsSolveRsig <- function(aug, params, ev, times, tol = 1e-10) {
    .eta <- params[["ETA_1_"]]
    list(
      f = c(1, 2),
      a = matrix(c(0.1, 0.2), 2L, 1L),
      A = array(0.01 + 0.5 * .eta, c(2L, 1L, 1L)),
      R = c(0.3, 0.4),
      aR = matrix(c(0.05, 0.06), 2L, 1L),
      AR = array(c(0.02, 0.03) + 0.25 * .eta, c(2L, 1L, 1L)),
      Rsig = matrix(c(0.7, 0.8), 2L, 1L),
      RsigDir = array(c(0.4, 0.45) + 0.125 * .eta, c(2L, 1L, 1L)),
      Rsig2 = array(c(0.6, 0.65) + 0.0625 * .eta, c(2L, 1L, 1L))
    )
  }

  test_that("the per-subject solve rebuilds the sigma slots the caller dropped (Route A)", {
    local_mocked_bindings(.foceiAnalyticSolveFA = .unitsSolveRsig)
    withr::local_envvar(FOCEI_NO_RSIG = "")
    .solve <- function(withR, sigSel) {
      .foceiAnalyticWhyReset()
      .foceiAnalyticSolveSubjectFD3(
        list(dirs = "ETA_1_"),
        c(ETA_1_ = 0),
        data.frame(),
        c(0, 1),
        withR = withR,
        sigSel = sigSel
      )
    }
    .r <- .solve(TRUE, 1L)
    # one direction plus the dropped sigma
    expect_identical(dim(.r$a), c(2L, 2L))
    expect_identical(.r$a[, 2], c(0, 0))
    expect_identical(.r$aR, cbind(c(0.05, 0.06), c(0.7, 0.8)))
    expect_equal(.r$AR[1, , ], rbind(c(0.02, 0.4), c(0.4, 0.6)), tolerance = 1e-12)
    expect_equal(.r$Ath[1, 1, , ], rbind(c(0.5, 0), c(0, 0)), tolerance = 1e-8)
    expect_equal(.r$AthR[1, 1, , ], rbind(c(0.25, 0.125), c(0.125, 0.0625)), tolerance = 1e-8)
    # without withR the slots are rebuilt too, with no third-order R
    .r <- .solve(FALSE, 1L)
    expect_identical(dim(.r$a), c(2L, 2L))
    expect_equal(.r$AR[1, , ], rbind(c(0.02, 0.4), c(0.4, 0.6)), tolerance = 1e-12)
    expect_null(.r$AthR)
    # a caller that dropped no sigma (the add/prop and SD-IOV assemblies) gets its own
    # directions only, R slots included
    .r <- .solve(TRUE, NULL)
    expect_identical(dim(.r$a), c(2L, 1L))
    expect_identical(dim(.r$AR), c(2L, 1L, 1L))
    expect_equal(drop(.r$AthR[1, 1, , ]), 0.25, tolerance = 1e-8)
    expect_identical(.r$R, c(0.3, 0.4))
    # FOCEI_NO_RSIG turns Route A off
    withr::local_envvar(FOCEI_NO_RSIG = "1")
    expect_identical(dim(.solve(TRUE, 1L)$a), c(2L, 1L))
  })

  # the same solve with a second residual sigma whose slots differ from the first's
  .unitsSolveRsig2 <- function(aug, params, ev, times, tol = 1e-10) {
    .E <- .unitsSolveRsig(aug, params, ev, times, tol)
    .eta <- params[["ETA_1_"]]
    .E$Rsig <- cbind(.E$Rsig, c(0.9, 1.1))
    .E$RsigDir <- array(c(.E$RsigDir, c(0.2, 0.25) + 0.5 * .eta), c(2L, 1L, 2L))
    .r2 <- array(0, c(2L, 2L, 2L))
    .r2[, 1, 1] <- .E$Rsig2[, 1, 1]
    .r2[, 2, 2] <- c(0.3, 0.35) + 0.75 * .eta
    .r2[, 1, 2] <- .r2[, 2, 1] <- c(0.05, 0.06)
    .E$Rsig2 <- .r2
    .E
  }

  test_that("the per-subject solve rebuilds only the sigma the caller dropped", {
    local_mocked_bindings(.foceiAnalyticSolveFA = .unitsSolveRsig2)
    withr::local_envvar(FOCEI_NO_RSIG = "")
    .foceiAnalyticWhyReset()
    # the second of two sigmas: its columns, not the first's
    .r <- .foceiAnalyticSolveSubjectFD3(
      list(dirs = "ETA_1_"),
      c(ETA_1_ = 0),
      data.frame(),
      c(0, 1),
      withR = TRUE,
      sigSel = 2L
    )
    expect_identical(dim(.r$a), c(2L, 2L))
    expect_identical(.r$aR, cbind(c(0.05, 0.06), c(0.9, 1.1)))
    expect_equal(.r$AR[1, , ], rbind(c(0.02, 0.2), c(0.2, 0.3)), tolerance = 1e-12)
    expect_equal(.r$AthR[1, 1, , ], rbind(c(0.25, 0.5), c(0.5, 0.75)), tolerance = 1e-8)
    # both, in the caller's order
    .r <- .foceiAnalyticSolveSubjectFD3(
      list(dirs = "ETA_1_"),
      c(ETA_1_ = 0),
      data.frame(),
      c(0, 1),
      withR = TRUE,
      sigSel = c(2L, 1L)
    )
    expect_identical(dim(.r$a), c(2L, 3L))
    expect_equal(.r$AR[1, , ], rbind(c(0.02, 0.2, 0.4), c(0.2, 0.3, 0.05), c(0.4, 0.05, 0.6)), tolerance = 1e-12)
  })

  test_that(".foceiAnalyticFD3Alloc(), Store() and Richardson() keep each slot apart", {
    .t <- .foceiAnalyticFD3Alloc(2L, 2L, 1L, 1L, TRUE, TRUE)
    expect_identical(names(.t), c("Ath", "AthR", "RsigDirEta", "Rsig2Eta"))
    expect_identical(dim(.t$RsigDirEta), c(2L, 2L, 1L, 1L))
    expect_null(.foceiAnalyticFD3Alloc(2L, 2L, 1L, 0L, TRUE, FALSE)$RsigDirEta)
    expect_null(.foceiAnalyticFD3Alloc(2L, 2L, 1L, 0L, FALSE, FALSE)$AthR)
    .d <- list(
      A = array(1, c(2, 1, 1)),
      AR = array(2, c(2, 1, 1)),
      RsigDir = array(3, c(2, 1, 1)),
      Rsig2 = array(4, c(2, 1, 1))
    )
    .t <- .foceiAnalyticFD3Store(.t, 2L, .d)
    expect_identical(
      c(.t$Ath[, 2, , ], .t$AthR[, 2, , ], .t$RsigDirEta[, 2, , ], .t$Rsig2Eta[, 2, , ]),
      rep(c(1, 2, 3, 4), each = 2)
    )
    expect_identical(c(.t$Ath[, 1, , ], .t$AthR[, 1, , ]), rep(0, 4))
    # x^3 at 1: the central difference at h is 3 + h^2, and (4 d(h/2) - d(h)) / 3
    # removes the h^2 term
    .s <- function(h) list(Ep = list(list(A = (1 + h)^3)), Em = list(list(A = (1 - h)^3)), h = h)
    expect_equal(.foceiAnalyticRichardson(.s(0.1), .s(0.1), 1L, "A"), 3 + 0.1^2, tolerance = 1e-12)
    expect_equal(.foceiAnalyticRichardson(.s(0.1), .s(0.05), 1L, "A"), 3, tolerance = 1e-12)
  })

  test_that(".foceiAnalyticRouteA() counts Rsig columns only for a caller that dropped sigmas", {
    withr::local_envvar(FOCEI_NO_RSIG = "")
    .E <- list(Rsig = matrix(0, 2L, 3L))
    expect_identical(.foceiAnalyticRouteA(.E, 2L), 3L)
    expect_identical(.foceiAnalyticRouteA(.E, NULL), 0L)
    expect_identical(.foceiAnalyticRouteA(.E, integer(0)), 0L)
    expect_identical(.foceiAnalyticRouteA(list(), 2L), 0L)
    expect_identical(.foceiAnalyticRouteA(list(Rsig = matrix(0, 2L, 0L)), 2L), 0L)
    withr::local_envvar(FOCEI_NO_RSIG = "1")
    expect_identical(.foceiAnalyticRouteA(.E, 2L), 0L)
  })

  test_that(".foceiAnalyticFD3Unflat() inverts .foceiAnalyticFD3Flat()", {
    .E <- list(
      A = array(1:8, c(2L, 2L, 2L)),
      AR = array(9:16, c(2L, 2L, 2L)),
      RsigDir = array(17:20, c(2L, 2L, 1L)),
      Rsig2 = array(21:22, c(2L, 1L, 1L))
    )
    for (.withR in c(FALSE, TRUE)) {
      for (.rsig in c(FALSE, TRUE)) {
        .v <- .foceiAnalyticFD3Flat(.E, .withR, .rsig)
        .u <- .foceiAnalyticFD3Unflat(.v, 2L, 2L, 1L, .withR, .rsig)
        .slots <- c("A", if (.withR) "AR", if (.rsig) c("RsigDir", "Rsig2"))
        expect_identical(names(.u), .slots)
        expect_identical(.u, .E[.slots])
      }
    }
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
