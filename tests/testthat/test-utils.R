nmTest({
  test_that(".collectWarn", {
    expect_equal(
      .collectWarn(
        {
          warning("A")
          1
        },
        lst = TRUE
      ),
      list(1, warning = "A", error = NULL)
    )
    # warnings are not duplicated
    expect_equal(
      .collectWarn(
        {
          warning("A")
          warning("A")
          1
        },
        lst = TRUE
      ),
      list(1, warning = "A", error = NULL)
    )
    # Non-lst raises the warning
    expect_warning(
      check <- .collectWarn(
        {
          warning("A")
          warning("A")
          1
        }
      ),
      regexp = "A"
    )
    expect_equal(check, 1)
    # Errors are not captured by default
    expect_error(
      .collectWarn(
        {
          warning("A")
          stop("B")
          1
        }
      ),
      regexp = "B"
    )
    # Errors are not captured by default
    expect_error(
      check <- .collectWarn(
        {
          warning("A")
          stop("B")
          1
        },
        lst = TRUE,
        collectErr = TRUE
      ),
      regexp = NA
    )
    expect_equal(check, list(NULL, warning = "A", error = "B"))
    # Inner try() blocks still work -- when expr ultimately succeeds, the
    # captured "internal" errors are discarded so they cannot leak into
    # the caller's error list (regression: previously every error fired
    # inside expr was returned, breaking SAEM's covariance recovery).
    expect_equal(
      .collectWarn(
        {
          try(stop("internal"), silent = TRUE)
          42
        },
        lst = TRUE,
        collectErr = TRUE
      ),
      list(42, warning = NULL, error = NULL)
    )
    # When the expression fails after an on.exit handler raises a second
    # error (issue 607 pattern), every message observed along the error
    # chain is preserved.  on.exit must be inside a function frame for
    # this to trigger; that mirrors the real-world rxUiGet.foceiHdEta
    # case in R/focei.R.
    .foceiHdEtaLike <- function() {
      on.exit(stop("Aborted calculation"))
      stop("none of the predictions depend on 'ETA'")
    }
    expect_error(
      check <- .collectWarn(.foceiHdEtaLike(), lst = TRUE, collectErr = TRUE),
      regexp = NA
    )
    expect_true(is.null(check[[1]]))
    expect_true("Aborted calculation" %in% check$error)
    expect_true("none of the predictions depend on 'ETA'" %in% check$error)
  })

  test_that(".nmObsRows() counts the records etTrans() keeps as observations", {
    .d <- data.frame(
      TIME = c(0, 1, 2, 3, 4, 5, 6, 7),
      EVID = c(1, 0, 0, 0, 2, 0, NA, 0),
      DV = c(NA, 1, NA, 3, NA, 4, 5, 6),
      MDV = c(1, 0, 0, 1, 1, NA, 0, 2)
    )
    # a missing EVID is 0 and only MDV = 1 is a missing observation, as in etTrans()
    expect_identical(.nmObsRows(.d), c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
    .d$MDV <- NULL
    expect_identical(.nmObsRows(.d), c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE))
    expect_identical(.nmObsRows(.d[0, ]), logical(0))
    # a missing TIME is no record at all to etTrans()
    .d$TIME[2] <- NA
    expect_false(.nmObsRows(.d)[2])
    # with no DV, a nonzero AMT, negative too, keeps an EVID 0 or missing-EVID
    # record (etTrans() drops its dose)
    .d <- data.frame(TIME = 1:6, EVID = c(NA, NA, NA, 0, 0, 0), DV = NA, AMT = c(50, -50, 0, 50, -50, 0))
    expect_identical(.nmObsRows(.d), c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE))
    # DV and TIME are read as FOCEi reads them, so "." is a missing DV
    .d <- data.frame(TIME = c("0", "1", "2"), EVID = 0, DV = c(".", "3", NA))
    expect_identical(.nmObsRows(.d), c(FALSE, TRUE, FALSE))
    # the columns are matched exactly: DVID is not DV
    expect_error(.nmObsRows(data.frame(TIME = 0, EVID = 0, DVID = 1)), "needs DV, EVID and TIME columns")
    expect_error(.nmObsRows(data.frame(TIME = 0, DV = 1)), "needs DV, EVID and TIME columns")
    expect_error(.nmObsRows(data.frame(EVID = 0, DV = 1)), "needs DV, EVID and TIME columns")
  })

  test_that(".nmObsRows() agrees with etTrans() over every EVID, DV, MDV, AMT and TIME case", {
    skip_on_cran()
    .m <- rxode2::rxode2({
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
    })
    .g <- expand.grid(
      evid = c(0, 1, 2, 3, 4, 101, NA),
      dv = c(4, NA),
      mdv = c(-1, 0, 1, NA, 2), # -1: no MDV column
      amt = c(0, NA, 50, -50),
      time = c(1, NA)
    )
    .n <- 0L
    for (.k in seq_len(nrow(.g))) {
      .r <- .g[.k, ]
      # a dose, the record under test, and an observation at t = 2
      .d <- data.frame(
        ID = 1L,
        TIME = c(0, .r$time, 2),
        DV = c(NA, .r$dv, 3),
        AMT = c(100, .r$amt, 0),
        EVID = c(1, .r$evid, 0),
        CMT = 2L
      )
      if (!(.r$mdv %in% -1)) {
        .d$MDV <- c(1, .r$mdv, 0)
      }
      # etTrans() rejects a dose with a missing AMT, and EVID 1 or 4 with an MDV of 0
      # or missing; those combinations are skipped, as no data can carry them
      .et <- tryCatch(as.data.frame(suppressWarnings(rxode2::etTrans(.d, .m))), error = function(e) NULL)
      if (is.null(.et)) {
        next
      }
      .n <- .n + 1L
      .want <- sum(.et$EVID == 0) == 2L
      expect_identical(.nmObsRows(.d)[2], .want, label = paste(names(.r), .r, collapse = " "))
    }
    # 560 combinations, of which etTrans() rejects 54
    expect_identical(.n, 506L)
  })

  test_that(".nmAssertIds() and .nmDropNoObsSubjects() match the FOCEi data step", {
    .d <- data.frame(ID = c(1, 1, 2, 2, 3), TIME = c(0, 1, 0, 1, 0), EVID = c(1, 0, 1, 0, 1), DV = c(NA, 2, NA, NA, NA))
    expect_warning(.k <- .nmDropNoObsSubjects(.d), "IDs without observations dropped: 2 3", fixed = TRUE)
    expect_identical(.k, .d[1:2, ])
    # all subjects have observations: unchanged; none has: kept whole
    expect_silent(expect_identical(.nmDropNoObsSubjects(.d[1:2, ]), .d[1:2, ]))
    expect_silent(expect_identical(.nmDropNoObsSubjects(.d[3:5, ]), .d[3:5, ]))
    # an observation with a missing TIME does not count, as in etTrans()
    .t <- data.frame(ID = c(1, 1, 2), TIME = c(0, 1, NA), EVID = c(1, 0, 0), DV = c(NA, 3, 4))
    expect_warning(.k <- .nmDropNoObsSubjects(.t), "IDs without observations dropped: 2", fixed = TRUE)
    expect_identical(.k, .t[1:2, ])
    .d$ID[2] <- NA
    expect_error(.nmDropNoObsSubjects(.d), "ID is missing on 1 record(s)", fixed = TRUE)
    expect_error(.nmSubjectObs(.d, 1), "ID is missing on 1 record(s)", fixed = TRUE)
    # the columns are matched exactly: IDENT is not ID
    .i <- data.frame(IDENT = 1, TIME = 0, EVID = 0, DV = 1)
    expect_error(.nmDropNoObsSubjects(.i), "the data needs an ID column", fixed = TRUE)
    expect_error(.nmSubjectObs(.i, 1), "the data needs an ID column", fixed = TRUE)
    expect_error(.nmSubjectObs(data.frame(ID = 1, EVID = 0, DV = 1), 1), "the data needs a TIME column", fixed = TRUE)
  })

  test_that(".nmDefaultEvid() makes a nonzero AMT a dose", {
    expect_identical(.nmDefaultEvid(data.frame(AMT = c(0, NA, 100, -1))), c(0L, 0L, 1L, 1L))
    expect_identical(.nmDefaultEvid(data.frame(DV = 1:2)), c(0L, 0L))
  })

  test_that(".nmSubjectObs() gives each subject's records and observations, in order", {
    .d <- data.frame(
      ID = c(2, 2, 2, 1, 1, 1),
      TIME = c(0, 1, 1, 0, 2, 3),
      EVID = c(1, 0, 0, 1, 0, 0),
      DV = c(NA, 5, NA, NA, 7, 8),
      CENS = c(0, 1, 0, 0, 0, 0)
    )
    .s <- .nmSubjectObs(.d, c(1, 2, 3))
    expect_length(.s, 3L)
    expect_null(.s[[3]])
    expect_identical(.s[[1]]$ev, .d[4:6, ])
    expect_identical(.s[[1]]$times, c(2, 3))
    expect_identical(.s[[1]]$y, c(7, 8))
    expect_identical(.s[[1]]$n, 2L)
    expect_identical(.s[[2]]$times, 1)
    expect_identical(.s[[2]]$y, 5)
    expect_identical(.s[[2]]$cens, 1L)
    expect_identical(.s[[2]]$limit, NA_real_)
    .d$CENS <- NULL
    expect_identical(.nmSubjectObs(.d, 1)[[1]]$cens, c(0L, 0L))
    # a subject's records come back in time order, which is the order a solve returns
    .u <- .d[c(6, 4, 5), ]
    .s <- .nmSubjectObs(.u, 1)[[1]]
    expect_identical(.s$ev, .d[4:6, ])
    expect_identical(.s$times, c(2, 3))
    expect_identical(.s$y, c(7, 8))
  })
})
