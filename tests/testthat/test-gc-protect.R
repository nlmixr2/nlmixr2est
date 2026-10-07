# Rcpp's create() and an Rcpp::Function call evaluate their arguments before
# they allocate (the result, or the call's pairlist), so a positional argument
# that is a freshly allocated, unprotected SEXP -- wrap(), dfCbindList() or
# getDfSubsetVars() -- can be collected by that allocation.  (A named SEXP,
# `_["x"] = wrap(y)`, is preserved by Rcpp; the scan flags it too, to keep one
# rule.)  Pass the C++ object itself, which is wrapped into the protected
# result, or hold the SEXP in an Rcpp object first.

.gcFreshSexp <- "\\b(wrap|dfCbindList|getDfSubsetVars)\\s*\\("

# The argument list of the balanced call whose opening parenthesis is just before
# character `open` of `txt`
.gcArgSpan <- function(open, txt) {
  .chars <- strsplit(substr(txt, open, open + 20000L), "", fixed = TRUE)[[1]]
  .depth <- 1L + cumsum((.chars == "(") - (.chars == ")"))
  paste(.chars[seq_len(which(.depth == 0L)[1] - 1L)], collapse = "")
}

# The `file:line` of every call in `txt` matching `callRe` with such an argument
.gcCallHits <- function(file, txt, callRe) {
  .starts <- gregexpr(callRe, txt, perl = TRUE)[[1]]
  if (.starts[1] < 0) {
    return(character(0))
  }
  .spans <- vapply(.starts + attr(.starts, "match.length"), .gcArgSpan, character(1), txt = txt)
  .hit <- grepl(.gcFreshSexp, .spans, perl = TRUE)
  if (!any(.hit)) {
    return(character(0))
  }
  .nl <- gregexpr("\n", txt, fixed = TRUE)[[1]]
  paste0(basename(file), ":", findInterval(.starts[.hit], .nl) + 1L)
}

# The create() calls, and the calls of the file's Rcpp::Function variables, in
# `file` with such an argument
.gcCreateHits <- function(file) {
  .txt <- paste(readLines(file, warn = FALSE), collapse = "\n")
  .create <- .gcCallHits(
    file,
    .txt,
    "\\b(List|DataFrame|CharacterVector|NumericVector|IntegerVector|LogicalVector)::create\\s*\\("
  )
  .fn <- regmatches(.txt, gregexpr("\\bFunction\\s+\\w+\\s*[=(]", .txt, perl = TRUE))[[1]]
  .fn <- unique(sub("^Function\\s+(\\w+).*", "\\1", .fn, perl = TRUE))
  if (length(.fn) == 0L) {
    return(.create)
  }
  c(
    .create,
    .gcCallHits(file, .txt, paste0("(?<![\\w.>])(", paste(.fn, collapse = "|"), ")\\s*\\("))
  )
}

test_that("no create() or Function call takes a freshly allocated SEXP as an argument", {
  .src <- testthat::test_path("..", "..", "src")
  skip_if_not(dir.exists(.src), "package sources not available")
  .files <- list.files(.src, pattern = "\\.(cpp|h)$", full.names = TRUE)
  expect_gt(length(.files), 10L)
  expect_identical(unlist(lapply(.files, .gcCreateHits)), character(0))
})

test_that("the scan finds a wrap() argument of create() or a Function, also across lines", {
  .f <- withr::local_tempfile(fileext = ".cpp")
  writeLines(
    c(
      "List a() { return List::create(_[\"x\"] = x); }",
      "List b() {",
      "  return List::create(_[\"x\"] = x,",
      "                      _[\"y\"] = wrap(y));",
      "}",
      "void c() {",
      "  Function covInstall = ns[\".f\"];",
      "  covInstall(e, xR);",
      "  covInstall(e,",
      "             wrap(y));",
      "}"
    ),
    .f
  )
  expect_identical(.gcCreateHits(.f), paste0(basename(.f), c(":3", ":9")))
})
