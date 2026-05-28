run_nlmixr2_runtime_benchmarks <- function(threads = NULL,
                                           file = NULL,
                                           calcTables = TRUE,
                                           quiet = TRUE) {
  nlmixr2est:::.nlmixr2BenchmarkRun(
    cases = nlmixr2est:::.nlmixr2BenchmarkDefaultCases(calcTables = calcTables),
    threads = threads,
    file = file,
    quiet = quiet
  )
}
