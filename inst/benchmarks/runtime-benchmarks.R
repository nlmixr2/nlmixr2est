run_nlmixr2_runtime_benchmarks <- function(threads = NULL,
                                           file = NULL,
                                           summary_file = NULL,
                                           calcTables = TRUE,
                                           quiet = TRUE,
                                           reps = 1L,
                                           warmup = 0L,
                                           summarize = FALSE) {
  .cases <- nlmixr2est:::.nlmixr2BenchmarkDefaultCases(calcTables = calcTables)
  .needsReplicates <- as.integer(reps) > 1L || as.integer(warmup) > 0L ||
    !is.null(summary_file) || isTRUE(summarize)
  if (!.needsReplicates) {
    return(nlmixr2est:::.nlmixr2BenchmarkRun(
      cases = .cases,
      threads = threads,
      file = file,
      quiet = quiet
    ))
  }
  .raw <- nlmixr2est:::.nlmixr2BenchmarkReplicate(
    cases = .cases,
    threads = threads,
    reps = reps,
    warmup = warmup,
    file = file,
    quiet = quiet
  )
  .summary <- nlmixr2est:::.nlmixr2BenchmarkSummarize(.raw)
  if (!is.null(summary_file)) {
    utils::write.csv(.summary, file = summary_file, row.names = FALSE)
  }
  if (isTRUE(summarize)) {
    return(.summary)
  }
  .raw
}
