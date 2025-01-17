library(testthat)
library(data.table)
library(rxode2)
library(nlmixr2est)
library(testthat)
if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
  # when testing CRAN, only use one thread
  setRxThreads(1L)
  setDTthreads(1L)
}

test_check("nlmixr2est")
## test_check("nlmixr2est", stop_on_failure = FALSE, wrap=TRUE,
##            reporter = testthat::LocationReporter)
