# Superassignment (`<<-`, `->>`) writes to a variable that is not named at the
# assignment site; keep an explicit environment (`acc$x <- ...`) instead.
# The parse data is inspected so strings and comments are not counted.
test_that("no R source file uses superassignment", {
  .rDir <- testthat::test_path("..", "..", "R")
  # an installed package (covr, R CMD check) has an R/ directory without sources
  .files <- list.files(.rDir, pattern = "\\.[Rr]$", full.names = TRUE)
  skip_if(length(.files) == 0L, "package source tree not available")
  .hits <- character(0)
  for (.f in .files) {
    .pd <- utils::getParseData(parse(.f, keep.source = TRUE))
    .bad <- .pd[
      (.pd$token == "LEFT_ASSIGN" & .pd$text == "<<-") |
        (.pd$token == "RIGHT_ASSIGN" & .pd$text == "->>"),
    ]
    if (nrow(.bad) > 0L) {
      .hits <- c(.hits, paste0(basename(.f), ":", .bad$line1))
    }
  }
  expect_identical(.hits, character(0))
})

test_that("the superassignment scan sees both operators", {
  .pd <- utils::getParseData(parse(
    text = c("f <- function() { a <<- 1; 2 ->> b; x <- '<<-' } # <<-"),
    keep.source = TRUE
  ))
  expect_identical(sum(.pd$token == "LEFT_ASSIGN" & .pd$text == "<<-"), 1L)
  expect_identical(sum(.pd$token == "RIGHT_ASSIGN" & .pd$text == "->>"), 1L)
})
