# Skip on CRAN only. CI's Windows and macOS runners take the CRAN subset
# (NOT_CRAN unset) and set CI=true, so a test inside
# nmTest(skipOnCran = FALSE) that calls this runs on every CI runner.
.skipOnCranOnly <- function() {
  testthat::skip_if(
    !identical(Sys.getenv("NOT_CRAN"), "true") && !isTRUE(as.logical(Sys.getenv("CI", "false"))),
    "on CRAN"
  )
}
