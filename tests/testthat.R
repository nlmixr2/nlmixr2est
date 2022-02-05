library(testthat)
library(rxode2)
library(nlmixr2)
verbose_minimization <- FALSE

test_check("nlmixr", stop_on_failure = FALSE, wrap=TRUE,
           reporter = testthat::LocationReporter)


