library(testthat)
library(data.table)
library(rxode2)
library(nlmixr2est)
library(testthat)
if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
  # when testing CRAN, only use two thread
  setRxThreads(2L)
  setDTthreads(2L)
  Sys.setenv(OMP_NUM_THREADS = "2")
  Sys.setenv(MKL_NUM_THREADS = "2")
  if (identical(Sys.info()["sysname"], "Darwin")) {
    rxode2::rxUnloadAll(set=FALSE)
  }
}

## test_check("nlmixr2est")
test_check("nlmixr2est", stop_on_failure = FALSE,
           reporter = testthat::LocationReporter)
