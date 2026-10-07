# Rcpp's create() evaluates its arguments before it allocates the result, so an
# argument that is a freshly allocated, unprotected SEXP -- wrap(), dfCbindList()
# or getDfSubsetVars() -- can be collected by that allocation.  Passing the C++
# object itself lets create() wrap it once the result is protected; anything else
# is held in an Rcpp object first.  This scans the sources for the pattern.

# The argument list of the balanced call whose opening parenthesis is just before
# character `open` of `txt`
.gcArgSpan <- function(open, txt) {
  .chars <- strsplit(substr(txt, open, open + 20000L), "", fixed = TRUE)[[1]]
  .depth <- 1L + cumsum((.chars == "(") - (.chars == ")"))
  paste(.chars[seq_len(which(.depth == 0L)[1] - 1L)], collapse = "")
}

# The `file:line` of every create() call in `file` with such an argument
.gcCreateHits <- function(file) {
  .txt <- paste(readLines(file, warn = FALSE), collapse = "\n")
  .starts <- gregexpr(
    "\\b(List|DataFrame|CharacterVector|NumericVector|IntegerVector|LogicalVector)::create\\s*\\(",
    .txt,
    perl = TRUE
  )[[1]]
  if (.starts[1] < 0) {
    return(character(0))
  }
  .spans <- vapply(.starts + attr(.starts, "match.length"), .gcArgSpan, character(1), txt = .txt)
  .hit <- grepl("\\b(wrap|dfCbindList|getDfSubsetVars)\\s*\\(", .spans, perl = TRUE)
  if (!any(.hit)) {
    return(character(0))
  }
  .nl <- gregexpr("\n", .txt, fixed = TRUE)[[1]]
  paste0(basename(file), ":", findInterval(.starts[.hit], .nl) + 1L)
}

test_that("no create() call takes a freshly allocated SEXP as an argument", {
  .src <- testthat::test_path("..", "..", "src")
  skip_if_not(dir.exists(.src), "package sources not available")
  .files <- list.files(.src, pattern = "\\.(cpp|h)$", full.names = TRUE)
  expect_gt(length(.files), 10L)
  expect_identical(unlist(lapply(.files, .gcCreateHits)), character(0))
})

test_that("the create() scan finds a wrap() argument, also across lines", {
  .f <- withr::local_tempfile(fileext = ".cpp")
  writeLines(
    c(
      "List a() { return List::create(_[\"x\"] = x); }",
      "List b() {",
      "  return List::create(_[\"x\"] = x,",
      "                      _[\"y\"] = wrap(y));",
      "}"
    ),
    .f
  )
  expect_identical(.gcCreateHits(.f), paste0(basename(.f), ":3"))
})
