.nlmixr <- function(...) {
  suppressWarnings(suppressMessages(nlmixr(...)))
}

# Drop-in replacement for set.seed() inside tests.  Sets the R RNG seed exactly
# like set.seed(), but also snapshots the R and rxode2 RNG state and restores
# both when the enclosing test (or calling frame) exits.  This keeps a test's
# data-sim seed from leaking forward and silently changing the data -- and thus
# the results -- of later tests that do not set their own seed.
.testSeed <- function(seed, envir = parent.frame()) {
  .oldR <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  .oldRx <- rxode2::rxGetSeed()
  withr::defer({
    if (is.null(.oldR)) {
      suppressWarnings(rm(".Random.seed", envir = .GlobalEnv))
    } else {
      assign(".Random.seed", .oldR, envir = .GlobalEnv)
    }
    rxode2::rxSetSeed(.oldRx)
  }, envir = envir)
  set.seed(seed)
  invisible()
}
# Use these for faster estimation (when numeric results don't matter but having
# an estimation result matters). Leave calcTables = TRUE since that is often
# needed for the tests.
saemControlFast <- saemControl(print = 0, nBurn = 1, nEm = 1, nmc = 1, nu = c(1, 1, 1))
foceiControlFast <- foceiControl(print = 0, maxInnerIterations = 1, maxOuterIterations = 1, eval.max = 1)
foceControlFast <- foceControl(print = 0, maxInnerIterations = 1, maxOuterIterations = 1, eval.max = 1)
