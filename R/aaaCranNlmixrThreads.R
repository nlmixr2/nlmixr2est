#' This function is to set the number of threads to 2
#'
#' In general it is a CRAN requirement that packages not use more than
#' 2 threads.  This function is to set the number of threads to 2 for
#' CRAN testing.  It is not intended for general use.
#'
#' When testing with devtools::test() or testthat::test_package(), the
#' NOT_CRAN environment variable is set to "true", so the number of
#' threads will not be limited to 2.
#'
#' @return nothing, called for side effect of setting the number of
#'   threads to 2 for CRAN testing
#'
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' # Set the number of threads to 2 for CRAN testing
#' aaaCranNlmixrThreads()
#'
aaaCranNlmixrThreads <- function() {
  if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
    # when testing CRAN, only use two thread
    rxode2::setRxThreads(2L)
    data.table::setDTthreads(2L)
    Sys.setenv(OMP_NUM_THREADS = "2")
    Sys.setenv(MKL_NUM_THREADS = "2")
    if (identical(Sys.info()["sysname"], "Darwin")) {
      rxode2::rxUnloadAll(set=FALSE)
    }
  }
}
