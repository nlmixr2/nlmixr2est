library(testthat)
library(data.table)
library(rxode2)
library(nlmixr2est)
library(testthat)
if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
  # when testing CRAN, only use one thread
  #setRxThreads(1L)
  #setDTthreads(1L)
  #rxode2::rxUnloadAll(FALSE)  # don't unload any models (seems to affect ASAN checks)
}

#rxCreateCache() ## Required for m1-san rhub

test_check("nlmixr2est")
## test_check("nlmixr2est", stop_on_failure = TRUE,
##            reporter = testthat::LocationReporter)
