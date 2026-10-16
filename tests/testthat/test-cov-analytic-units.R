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

  test_that("every exit of the entry points and assemblers is a named decline", {
    .other <- function(fn) {
      .k <- unname(.covAnalyticReturnKinds(fn))
      sort(.k[.k != ".foceiAnalyticFallback"])
    }
    expect_identical(.other(".foceiCovAnalyticCalc"), character(0))
    # the covType opt-out, which the C++ hook never reaches
    expect_identical(.other(".foceiCalcRanalyticAssemble"), "NULL")
    # the floored-R signal (twice) and the quadrature-node callback's NULL, which sets
    # the floored flag its caller turns into that signal
    expect_identical(.other(".foceiAnalyticAssembleR"), c("NULL", "structure", "structure"))
    # the FOCE branch's result
    expect_identical(.other(".foceiAnalyticAssembleRFR"), "R")
  })

  test_that("every helper of the decline chain that can give up is a leaf or exempt for a reason", {
    .reach <- .covAnalyticDeclineClosure()
    .leaves <- unique(unlist(.covAnalyticDeclineLeaves))
    expect_true(all(.leaves %in% .reach), label = paste(setdiff(.leaves, .reach), collapse = ", "))
    .gives <- .reach[vapply(.reach, .covAnalyticGivesUp, logical(1))]
    # a new helper that can give up must be made a leaf or exempted; a stale
    # exemption fails too
    expect_setequal(setdiff(.gives, .leaves), names(.covAnalyticDeclineExempt))
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
})
