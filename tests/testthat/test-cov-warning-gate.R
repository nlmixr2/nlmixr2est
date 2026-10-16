# A warning() that reports on a covariance method goes through .covWarn(), which
# clips it to the width of a $runInfo line.  The package's own functions are
# walked, so a warning added later is judged without being listed here.

# a warning whose text is about a covariance method has one of these in it
.covGatePattern <- paste(
  "covariance",
  "\\bFIM\\b",
  "sandwich",
  "\\b[RS] matrix",
  "\\$cov[A-Z]",
  "\\bcov(Full|Method)\\b",
  "\\bSE\\b",
  "\\banalytic err[0-9]",
  sep = "|"
)

.covGateWarners <- c("warning", "warningCondition", "simpleWarning", "warn", "cli_warn")

# The functions of `fns` that raise a covariance warning without .covWarn().  A
# warning() with no text of its own passes a condition on, so the text of its
# function stands for it.
.covGateBypasses <- function(fns, pattern = .covGatePattern, funnel = ".covWarn") {
  .hit <- character()
  for (.n in setdiff(names(fns), funnel)) {
    for (.call in .covGateFnCalls(fns[[.n]], .covGateWarners)) {
      .own <- .covGateStrings(.call)
      .txt <- if (length(.own) > 0L) .own else .covGateFnStrings(fns[[.n]])
      if (any(grepl(pattern, .txt))) {
        .hit <- c(.hit, .n)
      }
    }
  }
  unique(.hit)
}

# The functions of `fns` with a .covWarn() text whose fixed part alone is wider
# than a $runInfo line, labels and counts left out.
.covGateTooWide <- function(fns, width = .covWarnWidth) {
  .hit <- character()
  for (.n in names(fns)) {
    for (.call in .covGateFnCalls(fns[[.n]], ".covWarn")) {
      .fixed <- gsub("%%", "%", gsub("%[-+ 0#]*[0-9]*(\\.[0-9]+)?[sdif]", "", .covGateStrings(.call)))
      if (any(cli::ansi_nchar(.fixed, type = "width") > width)) {
        .hit <- c(.hit, .n)
      }
    }
  }
  unique(.hit)
}

# what a covariance warning is built or clipped with
.covGateTextFns <- c(".covWarn", ".covClip", ".covNotUsableMsg")

# The functions of the environment `ns` that call or pass by name one of `heads`
# (directly, in an argument default, as a symbol or as a string: do.call(".covWarn",
# ...)), the funnel apart.
.covGateCallers <- function(ns, heads = .covGateTextFns, funnel = ".covWarn") {
  .fns <- Filter(is.function, mget(ls(ns, all.names = TRUE), envir = ns))
  .has <- vapply(
    names(.fns),
    function(n) length(.covGateFnCalls(.fns[[n]], heads)) > 0L || any(heads %in% .covAnalyticNamedFns(n, ns)),
    logical(1)
  )
  setdiff(names(.fns)[.has], funnel)
}

# Callers a parse cannot see, C++ that warns with text built in R, and the R
# function that builds it.  foceiCalcRAnalytic (src/inner.cpp) warns with
# `.analyticWarnNoFallback`, which .foceiCalcRanalytic builds with
# .covNotUsableMsg (the C++ calls .covNotUsableMsg itself only when it is
# missing); its width is driven through that builder.
.covGateNative <- c(foceiCalcRAnalytic = ".foceiCalcRanalytic")

# The package functions the tests in `dir` say they drive: every string in the
# first argument of each .covSite() call.
.covGateDriven <- function(dir) {
  .files <- list.files(dir, pattern = "^test-.*[.]R$", full.names = TRUE)
  .names <- character()
  for (.f in .files) {
    for (.e in parse(.f, keep.source = FALSE)) {
      for (.call in .covGateWalk(.e, ".covSite")) {
        if (length(.call) >= 2L) {
          .names <- c(.names, .covGateStrings(.call[[2L]]))
        }
      }
    }
  }
  list(files = .files, names = unique(.names))
}

test_that("the gate finds a covariance warning that bypasses .covWarn()", {
  .bad <- list(
    plain = function() warning("covariance step failed"),
    sandwich = function() warning("since sandwich matrix is corrected", call. = FALSE),
    fim = function(x) warning(sprintf("linearized FIM failed; using %s", x)),
    rMatrix = function() warning(paste0("R matrix ", "is not finite")),
    sMatrix = function() warning("S matrix had problems"),
    covR = function() warning("compare to $covR"),
    covFull = function() warning("covFull variance block dropped"),
    covMethod = function() warning("switching to covMethod='r'"),
    se = function() warning("no SE for a theta"),
    cond = function() warningCondition("covariance step failed"),
    qualified = function() rlang::warn("covariance step failed"),
    nested = function() lapply(1, function(i) if (i) warning("covariance step failed")),
    default = function(a = warning("covariance step failed")) a,
    # no text of its own: the function's other text says what it is about
    passOn = function(w) {
      .why <- "covariance step failed"
      warning(w, call. = FALSE)
    }
  )
  expect_identical(.covGateBypasses(.bad), names(.bad))
})

test_that("the gate lets a covariance warning through .covWarn() and other warnings alone", {
  .ok <- list(
    viaFunnel = function() .covWarn("covariance step failed"),
    omega = function() warning("omega not positive definite; tables use its nearest PD matrix"),
    covariate = function() warning("covariate selection pinned to model-specified covariates"),
    covariates = function() warning("time-varying covariate(s) were excluded"),
    other = function() warning("outer bobyqa search stalled"),
    passOn = function(w) {
      .why <- "omega block zero"
      warning(w, call. = FALSE)
    },
    none = function() 1
  )
  expect_identical(.covGateBypasses(.ok), character())
  # the funnel is the one function allowed a warning() about a covariance
  .funnel <- list(.covWarn = function(w) {
    .why <- "covariance step failed"
    warning(w, call. = FALSE)
  })
  expect_identical(.covGateBypasses(.funnel), character())
  expect_identical(.covGateBypasses(.funnel, funnel = "other"), ".covWarn")
})

test_that("the gate lists the callers of the text functions and the ones a test drives", {
  .ns <- list2env(
    list(
      .covWarn = function(w) .covClip(w),
      .covClip = function(w, width = 74L) w,
      .covNotUsableMsg = function(request, why, tail) why,
      warns = function() .covWarn("x"),
      clips = function(w) .covClip(w, 54L),
      msg = function() .covNotUsableMsg("a", "b", "c"),
      nested = function() lapply(1, function(i) .covWarn("x")),
      default = function(a = .covClip("x")) a,
      byName = function() do.call(".covWarn", list("x")),
      passed = function(x) lapply(x, .covWarn),
      none = function() warning("x"),
      other = function() "a string"
    ),
    parent = emptyenv()
  )
  expect_setequal(.covGateCallers(.ns), c("warns", "clips", "msg", "nested", "default", "byName", "passed"))
  expect_setequal(
    .covGateCallers(.ns, funnel = "warns"),
    c(".covWarn", "clips", "msg", "nested", "default", "byName", "passed")
  )
  expect_setequal(.covGateCallers(.ns, heads = ".covNotUsableMsg"), "msg")
  .dir <- withr::local_tempdir()
  writeLines(
    c(
      "test_that('a', { .covSite('first', f()) })",
      "nmTest({ test_that('b', { x <- .covSite(\"second\", g()) }) })",
      "test_that('c', { .covSite(name, h()) })",
      "test_that('d', { .covSite(c('first', 'fourth'), f()) })"
    ),
    file.path(.dir, "test-one.R")
  )
  writeLines("test_that('e', { .covSite('third', f()) })", file.path(.dir, "test-two.R"))
  writeLines(".covSite('hidden', f())", file.path(.dir, "other.R"))
  .d <- .covGateDriven(.dir)
  expect_identical(basename(.d$files), c("test-one.R", "test-two.R"))
  expect_identical(.d$names, c("first", "second", "fourth", "third"))
})

test_that("every function that builds or clips a covariance warning is driven by a .covSite() test", {
  .d <- .covGateDriven(testthat::test_path())
  expect_gt(length(.d$files), 10L)
  .callers <- .covGateCallers(asNamespace("nlmixr2est"))
  # a caller no test drives, and a name no function answers to
  expect_identical(setdiff(.callers, .d$names), character())
  expect_identical(setdiff(.d$names, .callers), character())
  # the C++ callers get their text from an R function that is driven
  expect_true(all(.covGateNative %in% .callers))
  expect_true(all(.covGateNative %in% .d$names))
})

test_that("no warning in the package reports on a covariance method except through .covWarn()", {
  .fns <- .covGateFunctions()
  # the funnel itself is the one warning() call it makes, and the package uses it
  expect_length(.covGateFnCalls(.covWarn, "warning"), 1L)
  expect_true(any(vapply(.fns, function(f) length(.covGateFnCalls(f, ".covWarn")) > 0L, logical(1))))
  expect_identical(.covGateBypasses(.fns), character())
})

# a function whose body is `call`, with the text inside it as a literal
.covGateFn <- function(call) {
  eval(bquote(function(a = 1) .(call)))
}

test_that("the gate finds a .covWarn() text too wide before any label is added", {
  .x <- function(n) strrep("x", n)
  .wide <- list(
    long = .covGateFn(call(".covWarn", .x(75))),
    format = .covGateFn(call(".covWarn", call("sprintf", paste0("\"%s\" ", .x(72), " %d"), quote(a), 1L))),
    branch = .covGateFn(call(".covWarn", call("if", quote(a), "short", .x(80)))),
    nested = .covGateFn(call("lapply", 1, call("function", as.pairlist(alist(i = )), call(".covWarn", .x(75)))))
  )
  expect_identical(.covGateTooWide(.wide), names(.wide))
  .fit <- list(
    exact = .covGateFn(call(".covWarn", .x(74))),
    format = .covGateFn(call(".covWarn", call("sprintf", paste0("\"%s\" ", .x(70), " %d"), quote(a), 1L))),
    percent = .covGateFn(call(".covWarn", call("sprintf", paste0("100%% ", .x(69))))),
    dynamic = function(w) .covWarn(w)
  )
  expect_identical(.covGateTooWide(.fit), character())
})

test_that("every .covWarn() text in the package fits a $runInfo line before any label is added", {
  expect_identical(.covGateTooWide(.covGateFunctions()), character())
})
