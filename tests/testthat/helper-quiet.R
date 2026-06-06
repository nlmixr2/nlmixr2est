.nlmixr <- function(...) {
  suppressWarnings(suppressMessages(nlmixr(...)))
}
# Use these for faster estimation (when numeric results don't matter but having
# an estimation result matters). Leave calcTables = TRUE since that is often
# needed for the tests.
saemControlFast <- saemControl(print = 0, nBurn = 1, nEm = 1, nmc = 1, nu = c(1, 1, 1))
foceiControlFast <- foceiControl(print = 0, maxInnerIterations = 1, maxOuterIterations = 1, eval.max = 1)
foceControlFast <- foceControl(print = 0, maxInnerIterations = 1, maxOuterIterations = 1, eval.max = 1)
