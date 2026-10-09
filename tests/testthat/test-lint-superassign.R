# Superassignment (`<<-`, `->>`) writes to a variable that is not named at the
# assignment site; keep an explicit environment (`acc$x <- ...`) instead.
# The parse data is inspected so strings and comments are not counted.
test_that("no R source file uses superassignment", {
  .rDir <- testthat::test_path("..", "..", "R")
  # an installed package (covr, R CMD check) has an R/ directory too, holding the
  # lazy-load database instead of sources
  skip_if_not(file.exists(file.path(.rDir, "nlmixr2.R")), "package source tree not available")
  .files <- list.files(.rDir, pattern = "\\.[Rr]$", full.names = TRUE)
  expect_gt(length(.files), 0L)
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
    text = "f <- function() { a <<- 1; 2 ->> b; x <- '<<-' } # <<-",
    keep.source = TRUE
  ))
  expect_identical(sum(.pd$token == "LEFT_ASSIGN" & .pd$text == "<<-"), 1L)
  expect_identical(sum(.pd$token == "RIGHT_ASSIGN" & .pd$text == "->>"), 1L)
})
