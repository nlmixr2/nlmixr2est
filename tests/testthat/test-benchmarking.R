nmTest({
  test_that("benchmark timing normalization maps estimator stages", {
    .foceiTime <- data.frame(setup = 1, optimize = 2, covariance = 3, check.names = FALSE, row.names = "elapsed")
    .saemTime <- data.frame(
      preprocess = 0.5,
      setup = 1.5,
      saem = 2.5,
      covariance = 0.7,
      check.names = FALSE,
      row.names = "elapsed"
    )

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
        control = foceiControl(
          print = 0L,
          maxInnerIterations = 1L,
          maxOuterIterations = 1L,
          eval.max = 1L,
          calcTables = FALSE
        )
      )
    ))

    expect_identical(row.names(.fit$time), "elapsed")
    expect_true(all(c("preprocess", "setup", "optimize", "covariance", "compress") %in% names(.fit$time)))

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
          foceiControl(print = 0L, maxInnerIterations = 1L, maxOuterIterations = 1L, eval.max = 1L, calcTables = FALSE)
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

  test_that("the benchmark counts the observations etTrans() keeps", {
    # lower-case names, a missing DV, MDV = 1, a missing EVID with a DV, and a dose
    .d <- data.frame(
      id = 1,
      time = 0:5,
      evid = c(1, 0, 0, 0, NA, 2),
      dv = c(NA, 1, NA, 3, 4, NA),
      mdv = c(1, 0, 0, 1, 0, 1)
    )
    expect_identical(.nlmixr2BenchmarkNobs(.d), 2L)
    # no EVID: taken from AMT as the VAE prep does (a nonzero AMT is a dose)
    expect_identical(.nlmixr2BenchmarkNobs(.d[, c("id", "time", "dv")]), 3L)
    .a <- data.frame(id = 1, time = 0:3, dv = c(5, 1, 2, 3), amt = c(100, 0, 0, 0))
    expect_identical(.nlmixr2BenchmarkNobs(.a), 3L)
    # no DV or no TIME: nothing to count
    expect_identical(.nlmixr2BenchmarkNobs(.d[, c("id", "time", "evid")]), NA_integer_)
    expect_identical(.nlmixr2BenchmarkNobs(.d[, c("id", "evid", "dv")]), NA_integer_)
    expect_identical(.nlmixr2BenchmarkNobs(NULL), NA_integer_)
    expect_identical(.nlmixr2BenchmarkNobs(nlmixr2data::theo_sd), 132L)
  })

  test_that(".nlmixr2BenchmarkNobs() gives data without an EVID column FOCEi's EVID", {
    # no AMT column: every record EVID 0, so only a missing DV is not counted
    expect_identical(.nlmixr2BenchmarkNobs(data.frame(ID = 1, TIME = 0:2, DV = c(1, NA, 3))), 2L)
    # AMT and MDV: MDV is the EVID, so MDV 0 with an AMT is an observation
    .d <- data.frame(ID = 1, TIME = 0:2, DV = c(NA, 4, 5), AMT = c(100, 50, 0), MDV = c(1, 0, 0))
    expect_identical(.nlmixr2BenchmarkNobs(.d), 2L)
    # AMT only: a nonzero AMT is a dose
    .d$MDV <- NULL
    expect_identical(.nlmixr2BenchmarkNobs(.d), 1L)
  })
})
