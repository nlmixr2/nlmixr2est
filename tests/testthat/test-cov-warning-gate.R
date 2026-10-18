# A warning() that reports on a covariance method goes through .covWarn(), which
# clips it to the width of a $runInfo line.  The package's own functions and the
# C++ sources are walked, so a warning added later is judged without being listed
# here.

# A warning whose text is about a covariance method has one of these in it, in any
# case.  The words of the label-only style (a quoted method name, "kept",
# "installed", "not usable") are here because the covariance warnings no longer
# name a covariance.
.covGatePattern <- paste(
  "covariance",
  "\\bFIM\\b",
  "sandwich",
  "\\b[RS] matrix",
  "\\$cov[A-Z]",
  "\\bcov(Full|Method)\\b",
  "\\bSE\\b",
  "\\banalytic err[0-9]",
  "standard error",
  "information matrix",
  "sqrtm",
  "\\bkept\\b",
  "\\binstalled\\b",
  "not usable",
  "linFim",
  "\\bHa\\b",
  "^\"",
  sep = "|"
)

.covGateWarners <- c("warning", "warningCondition", "simpleWarning", "warn", "cli_warn")

# what builds the text of a covariance warning: a function that calls one of these
# and then warning() has built a warning that skips the clip
.covGateBuilders <- c(".covClip", ".covNotUsableMsg", ".covInsteadMsg")

# The functions of `fns` that raise a covariance warning without .covWarn().  A
# warning() with no text of its own passes a condition on, so the text of its
# function stands for it; a function that builds a text with a builder and also
# calls warning() is taken to pass that text on.
.covGateBypasses <- function(fns, pattern = .covGatePattern, funnel = ".covWarn", builders = .covGateBuilders) {
  .hit <- character()
  for (.n in setdiff(names(fns), funnel)) {
    .calls <- .covGateFnCalls(fns[[.n]], .covGateWarners)
    if (length(.calls) > 0L && length(.covGateFnCalls(fns[[.n]], builders)) > 0L) {
      .hit <- c(.hit, .n)
      next
    }
    for (.call in .calls) {
      .own <- .covGateStrings(.call)
      .txt <- if (length(.own) > 0L) .own else .covGateFnStrings(fns[[.n]])
      if (any(grepl(pattern, .txt, ignore.case = TRUE))) {
        .hit <- c(.hit, .n)
      }
    }
  }
  unique(.hit)
}

# a format string without its conversions, and %% as one %
.covGateStrip <- function(s) {
  gsub("%%", "%", gsub("%[-+ 0#]*[0-9]*(\\.[0-9]+)?[sdif]", "", s))
}

# The columns of text an expression always gives, whatever it is filled with: a
# string, the sum of the parts of paste0(), paste() and sprintf(), the longer
# branch of if and ifelse, the last expression of a block.
.covGateFixedWidth <- function(e) {
  if (is.character(e)) {
    return(sum(cli::ansi_nchar(.covGateStrip(e), type = "width")))
  }
  if (!is.call(e)) {
    return(0L)
  }
  .h <- .covGateHead(e)
  .args <- as.list(e)[-1L]
  .keep <- vapply(seq_along(.args), function(i) !identical(.args[[i]], quote(expr = )), logical(1))
  .args <- .args[.keep]
  .names <- names(.args)
  if (is.null(.names)) {
    .names <- rep("", length(.args))
  }
  if (.h == "{") {
    return(if (length(.args) > 0L) .covGateFixedWidth(.args[[length(.args)]]) else 0L)
  }
  if (.h %in% c("if", "ifelse")) {
    .branches <- .args[-1L]
    return(if (length(.branches) > 0L) max(vapply(.branches, .covGateFixedWidth, numeric(1))) else 0L)
  }
  if (.h %in% c("paste0", "paste", "sprintf", "gettextf")) {
    .parts <- .args[!(.names %in% c("sep", "collapse"))]
    .sum <- sum(vapply(.parts, .covGateFixedWidth, numeric(1)))
    if (.h == "paste") {
      .sep <- if ("sep" %in% .names && is.character(.args[[which(.names == "sep")[1L]]])) {
        .args[[which(.names == "sep")[1L]]]
      } else {
        " "
      }
      .sum <- .sum + max(length(.parts) - 1L, 0L) * nchar(.sep)
    }
    return(.sum)
  }
  0L
}

# The functions of `fns` with a .covWarn() text, or a text builder, whose fixed part
# alone is wider than a $runInfo line, labels and counts left out.
.covGateTooWide <- function(fns, width = .covWarnWidth, builders = ".covInsteadMsg") {
  .hit <- character()
  for (.n in names(fns)) {
    .wide <- vapply(
      .covGateFnCalls(fns[[.n]], ".covWarn"),
      function(.call) length(.call) >= 2L && .covGateFixedWidth(.call[[2L]]) > width,
      logical(1)
    )
    if (any(.wide) || (.n %in% builders && .covGateFixedWidth(body(fns[[.n]])) > width)) {
      .hit <- c(.hit, .n)
    }
  }
  unique(.hit)
}

# what a covariance warning is built or clipped with
.covGateTextFns <- c(".covWarn", ".covClip", ".covNotUsableMsg", ".covInsteadMsg")

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
# missing), and with `.covBoundaryWarn()`; their width is driven through those.
.covGateNative <- c(foceiCalcRAnalytic = ".foceiCalcRanalytic", foceiCalcCov = ".covBoundaryWarn")

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

# The text of a C++ call that opens at `start` (just after its parenthesis), up to
# the parenthesis that closes it; parentheses inside string literals do not count.
.covGateCppCall <- function(txt, start) {
  .chars <- strsplit(substr(txt, start, min(start + 2000L, nchar(txt))), "")[[1L]]
  .depth <- 1L
  .quote <- FALSE
  .i <- 1L
  while (.i <= length(.chars)) {
    .ch <- .chars[.i]
    if (.quote) {
      if (.ch == "\\") {
        .i <- .i + 1L
      } else if (.ch == "\"") {
        .quote <- FALSE
      }
    } else if (.ch == "\"") {
      .quote <- TRUE
    } else if (.ch == "(") {
      .depth <- .depth + 1L
    } else if (.ch == ")") {
      .depth <- .depth - 1L
      if (.depth == 0L) {
        return(paste(.chars[seq_len(.i - 1L)], collapse = ""))
      }
    }
    .i <- .i + 1L
  }
  stop("a C++ warning call that does not close", call. = FALSE)
}

# The warning() and Rf_warning() calls in the C++ sources of `dir` whose text
# mentions a covariance: the text the literals of the call give, its columns
# without conversions, and its lines.
.covGateCpp <- function(dir) {
  .files <- list.files(dir, pattern = "[.](cpp|c|h)$", full.names = TRUE)
  .out <- list()
  for (.f in .files) {
    .txt <- paste(readLines(.f, warn = FALSE), collapse = "\n")
    .at <- gregexpr("(?<![A-Za-z_0-9])(Rf_|Rcpp::)?warning\\s*\\(", .txt, perl = TRUE)[[1L]]
    if (.at[1L] < 0L) {
      next
    }
    for (.k in seq_along(.at)) {
      .call <- .covGateCppCall(.txt, .at[.k] + attr(.at, "match.length")[.k])
      .lit <- regmatches(.call, gregexpr("\"((?:[^\"\\\\]|\\\\.)*)\"", .call, perl = TRUE))[[1L]]
      .msg <- gsub("\\\\n", "\n", gsub("^\"|\"$", "", .lit))
      .msg <- .covGateStrip(paste(.msg, collapse = ""))
      if (grepl("covariance|getVarCov|\\bSE\\b|sandwich|\\bFIM\\b", .msg, ignore.case = TRUE)) {
        .out[[length(.out) + 1L]] <- data.frame(
          file = basename(.f),
          text = .msg,
          width = cli::ansi_nchar(.msg, type = "width"),
          lines = length(strsplit(.msg, "\n", fixed = TRUE)[[1L]])
        )
      }
    }
  }
  do.call(rbind, .out)
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
    },
    # the label-only style
    labelOnly = function(x) warning(sprintf("\"linFim\" %s; %s", x, x)),
    rejectCopy = function() warning("\"%s\" %s; none installed"),
    capital = function() warning("Covariance step failed"),
    standardErrors = function() warning("standard errors could not be computed"),
    informationMatrix = function() warning("the information matrix is singular"),
    repair = function() warning("repaired with sqrtm"),
    kept = function() warning("kept the previous one"),
    notUsable = function() warning("the request is not usable"),
    ha = function() warning("Ha is not finite"),
    # a text built by a builder and warned with a plain warning()
    notUsableThenWarning = function() {
      .m <- .covNotUsableMsg("a", "b", "c")
      warning(.m)
    },
    clipThenWarning = function(w) {
      .t <- .covClip(w)
      warning(.t, call. = FALSE)
    },
    insteadThenWarning = function(l, w) {
      .t <- .covInsteadMsg(l, w)
      warning(.t, call. = FALSE)
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
    hessian = function() warning("Outer Hessian unavailable; restarting gradient-only nlminb"),
    passOn = function(w) {
      .why <- "omega block zero"
      warning(w, call. = FALSE)
    },
    builderOnly = function() .covNotUsableMsg("a", "b", "c"),
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
      .covInsteadMsg = function(label, what) label,
      warns = function() .covWarn("x"),
      clips = function(w) .covClip(w, 54L),
      msg = function() .covNotUsableMsg("a", "b", "c"),
      instead = function() .covInsteadMsg("a", "b"),
      nested = function() lapply(1, function(i) .covWarn("x")),
      default = function(a = .covClip("x")) a,
      byName = function() do.call(".covWarn", list("x")),
      passed = function(x) lapply(x, .covWarn),
      none = function() warning("x"),
      other = function() "a string"
    ),
    parent = emptyenv()
  )
  .all <- c("warns", "clips", "msg", "instead", "nested", "default", "byName", "passed")
  expect_setequal(.covGateCallers(.ns), .all)
  expect_setequal(.covGateCallers(.ns, funnel = "warns"), c(".covWarn", setdiff(.all, "warns")))
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
  # walks every function of the package
  skip_on_cran()
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
  # walks every function of the package
  skip_on_cran()
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
  .p0 <- function(...) as.call(c(list(quote(paste0)), list(...)))
  .wide <- list(
    long = .covGateFn(call(".covWarn", .x(75))),
    format = .covGateFn(call(".covWarn", call("sprintf", paste0("\"%s\" ", .x(72), " %d"), quote(a), 1L))),
    branch = .covGateFn(call(".covWarn", call("if", quote(a), "short", .x(80)))),
    nested = .covGateFn(call("lapply", 1, call("function", as.pairlist(alist(i = )), call(".covWarn", .x(75))))),
    # the pieces of paste0() and paste() add up
    pasted = .covGateFn(call(".covWarn", .p0(.x(60), .x(60)))),
    pastedSep = .covGateFn(call(".covWarn", call("paste", .x(40), .x(34)))),
    pastedInSprintf = .covGateFn(call(".covWarn", call("sprintf", "%s done", .p0(.x(40), .x(30), quote(a))))),
    block = .covGateFn(call(".covWarn", call("{", quote(a), .x(80)))),
    ifelse = .covGateFn(call(".covWarn", call("ifelse", quote(a), .x(80), "short")))
  )
  expect_identical(.covGateTooWide(.wide), names(.wide))
  .fit <- list(
    exact = .covGateFn(call(".covWarn", .x(74))),
    format = .covGateFn(call(".covWarn", call("sprintf", paste0("\"%s\" ", .x(70), " %d"), quote(a), 1L))),
    percent = .covGateFn(call(".covWarn", call("sprintf", paste0("100%% ", .x(69))))),
    pasted = .covGateFn(call(".covWarn", .p0(.x(40), .x(34), quote(a)))),
    pastedSep = .covGateFn(call(".covWarn", call("paste", .x(40), .x(33)))),
    pastedCollapse = .covGateFn(call(".covWarn", call("paste0", .x(70), quote(a), collapse = .x(80)))),
    dynamic = function(w) .covWarn(w)
  )
  expect_identical(.covGateTooWide(.fit), character())
  # a text builder is measured where it is written
  .builders <- list(
    .covInsteadMsg = .covGateFn(call("sprintf", paste0("\"%s\" ", .x(75)), quote(a))),
    other = .covGateFn(call("sprintf", paste0("\"%s\" ", .x(75)), quote(a)))
  )
  expect_identical(.covGateTooWide(.builders), ".covInsteadMsg")
  expect_identical(.covGateTooWide(.builders["other"]), character())
})

test_that("every .covWarn() text in the package fits a $runInfo line before any label is added", {
  # walks every function of the package
  skip_on_cran()
  expect_identical(.covGateTooWide(.covGateFunctions()), character())
})

test_that("the gate reads the covariance warnings of C++ sources", {
  .dir <- withr::local_tempdir()
  writeLines(
    c(
      "void f() {",
      "  warning(_(\"covariance step failed\"));",
      "  Rf_warning(_(\"covariance failed; objective not evaluable at %d\"), 3);",
      "  warning(_(\"the SE (see $covMethod); a check\"));",
      "  // too wide, in two literals",
      "  warning(_(\"the covariance of the estimates is not available for \"",
      "            \"this fit because of a problem that has a long description\"));",
      "  Rcpp::warning(_(\"parameter estimate near boundary; covariance not calculated:\\n   \") + s);",
      "  warning(_(\"an optimizer message that is long but is about something else, not the thing\"));",
      "  stop(\"covariance failed\");",
      "}"
    ),
    file.path(.dir, "a.cpp")
  )
  writeLines("void g() { warning(_(\"gradient problems with covariance; see $scaleInfo\")); }", file.path(.dir, "b.h"))
  writeLines("warning(\"covariance, not in a source\")", file.path(.dir, "notes.txt"))
  .cpp <- .covGateCpp(.dir)
  expect_identical(.cpp$file, c(rep("a.cpp", 5L), "b.h"))
  expect_identical(
    .cpp$width[c(1L, 2L, 3L, 6L)],
    nchar(c(
      "covariance step failed",
      "covariance failed; objective not evaluable at ",
      "the SE (see $covMethod); a check",
      "gradient problems with covariance; see $scaleInfo"
    ))
  )
  expect_identical(.cpp$lines, c(1L, 1L, 1L, 1L, 2L, 1L))
  expect_identical(.cpp$text[5L], "parameter estimate near boundary; covariance not calculated:\n   ")
  # one literal too wide, and the one with two lines
  expect_identical(which(.cpp$width > .covWarnWidth | .cpp$lines > 1L), c(4L, 5L))
})

test_that("no C++ covariance warning is wider than a $runInfo line", {
  # the sources are not in an installed package
  .src <- testthat::test_path("..", "..", "src")
  skip_if_not(dir.exists(.src), "no C++ sources to read")
  .cpp <- .covGateCpp(.src)
  # the scan reads the warnings of covariance steps: they are all here
  expect_gt(nrow(.cpp), 8L)
  expect_true(all(
    c("covariance step failed", "analytic covariance: tolerances increased for some subjects") %in% .cpp$text
  ))
  expect_identical(.cpp[.cpp$width > .covWarnWidth | .cpp$lines > 1L, "text"], character())
})
