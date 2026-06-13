nmTest({
  test_that("benchmark timing normalization maps estimator stages", {
    .foceiTime <- data.frame(setup = 1, optimize = 2, covariance = 3,
                             check.names = FALSE, row.names = "elapsed")
    .saemTime <- data.frame(preprocess = 0.5, setup = 1.5, saem = 2.5, covariance = 0.7,
                            check.names = FALSE, row.names = "elapsed")

    .foceiStages <- .nlmixr2BenchmarkNormalizeTime(.foceiTime)
    .saemStages <- .nlmixr2BenchmarkNormalizeTime(.saemTime)

    expect_equal(.foceiStages$elapsed[.foceiStages$stage == "estimate"], 2)
    expect_equal(.saemStages$elapsed[.saemStages$stage == "estimate"], 2.5)
    expect_equal(.saemStages$elapsed[.saemStages$stage == "preprocess"], 0.5)
  })

  test_that("focei timing uses elapsed row names and benchmark extractor returns tidy rows", {
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(
        one.compartment,
        nlmixr2data::theo_sd,
        est = "focei",
        control = foceiControl(print = 0L, maxInnerIterations = 1L,
                               maxOuterIterations = 1L, eval.max = 1L,
                               calcTables = FALSE)
      )
    ))

    expect_identical(row.names(.fit$time), "elapsed")
    # `compress` only appears when its stage exceeds the 5e-5s timing floor
    # (.finalizeOverallTiming in R/timing.R).  With calcTables = FALSE it is ~0
    # and dropped -- reliably so on Windows' coarser clock -- so require only
    # the stages that always take measurable time.
    expect_true(all(c("preprocess", "setup", "optimize", "covariance") %in% names(.fit$time)))

    .bench <- .nlmixr2BenchmarkExtractFit(.fit, "fresh-focei")
    expect_true(all(c("case", "estimator", "stage", "raw_stage", "elapsed", "total_elapsed") %in% names(.bench)))
    expect_true(any(.bench$stage == "estimate"))
  })

  test_that("benchmark runner can execute a custom case and export csv", {
    .csv <- tempfile(fileext = ".csv")
    .cases <- list(
      fast_focei = list(
        est = "focei",
        model = one.compartment,
        data = nlmixr2data::theo_sd,
        control = function() {
          foceiControl(print = 0L, maxInnerIterations = 1L,
                       maxOuterIterations = 1L, eval.max = 1L,
                       calcTables = FALSE)
        }
      )
    )

    .bench <- .nlmixr2BenchmarkRun(cases = .cases, file = .csv, quiet = TRUE)

    expect_true(file.exists(.csv))
    expect_true(nrow(.bench) > 0)
    expect_true(any(.bench$stage == "estimate"))
  })

  test_that("benchmark summaries aggregate replicated runs", {
    .bench <- data.frame(
      case = c("x", "x", "x", "x"),
      estimator = c("focei", "focei", "focei", "focei"),
      method = c("FOCE", "FOCE", "FOCE", "FOCE"),
      subjects = c(12L, 12L, 12L, 12L),
      observations = c(132L, 132L, 132L, 132L),
      ntheta = c(4L, 4L, 4L, 4L),
      neta = c(3L, 3L, 3L, 3L),
      rx_threads = c(1L, 1L, 1L, 1L),
      omp_threads = c(NA_character_, NA_character_, NA_character_, NA_character_),
      total_elapsed = c(5, 7, 5, 7),
      stage = c("estimate", "estimate", "covariance", "covariance"),
      raw_stage = c("optimize", "optimize", "covariance", "covariance"),
      elapsed = c(4, 6, 1, 1),
      replicate = c(1L, 2L, 1L, 2L),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    .summary <- .nlmixr2BenchmarkSummarize(.bench)

    expect_equal(.summary$reps[.summary$stage == "estimate"], 2)
    expect_equal(.summary$elapsed_mean[.summary$stage == "estimate"], 5)
    expect_equal(.summary$total_elapsed_mean[.summary$stage == "covariance"], 6)
  })
})
