library(testthat)
library(rxode2)
library(nlmixr2est)
verbose_minimization <- FALSE
test_check("nlmixr2est")
## test_check("nlmixr2est", stop_on_failure = FALSE, wrap=TRUE,
##            reporter = testthat::LocationReporter)
