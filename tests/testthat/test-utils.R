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
      EVID = c(1, 0, 0, 0, 2, 0, NA, 0),
      DV = c(NA, 1, NA, 3, NA, 4, 5, 6),
      MDV = c(1, 0, 0, 1, 1, NA, 0, 2)
    )
    # a missing EVID is 0 and only MDV = 1 is a missing observation, as in etTrans()
    expect_identical(.nmObsRows(.d), c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
    .d$MDV <- NULL
    expect_identical(.nmObsRows(.d), c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE))
    expect_identical(.nmObsRows(.d[0, ]), logical(0))
    # etTrans() keeps a missing-EVID record with a positive AMT even with no DV
    .d <- data.frame(EVID = c(NA, NA), DV = c(NA, NA), AMT = c(50, 0))
    expect_identical(.nmObsRows(.d), c(TRUE, FALSE))
    # the columns are matched exactly: DVID is not DV
    expect_error(.nmObsRows(data.frame(EVID = 0, DVID = 1)), "needs a DV and an EVID column")
    expect_error(.nmObsRows(data.frame(DV = 1)), "needs a DV and an EVID column")
  })

  test_that(".nmObsRows() agrees with etTrans() over every EVID, DV, MDV and AMT case", {
    skip_on_cran()
    .m <- rxode2::rxode2({
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
    })
    .g <- expand.grid(
      evid = c(0, NA),
      dv = c(4, NA),
      mdv = c(-1, 0, 1, NA, 2), # -1: no MDV column
      amt = c(0, NA, 50)
    )
    .g <- .g[!(.g$evid %in% 0 & .g$amt %in% 50), ]
    for (.k in seq_len(nrow(.g))) {
      .r <- .g[.k, ]
      # a dose, the record under test at t = 1, and an observation
      .d <- data.frame(
        ID = 1L,
        TIME = c(0, 1, 2),
        DV = c(NA, .r$dv, 3),
        AMT = c(100, .r$amt, 0),
        EVID = c(1, .r$evid, 0),
        CMT = 2L
      )
      if (!(.r$mdv %in% -1)) {
        .d$MDV <- c(1, .r$mdv, 0)
      }
      .et <- as.data.frame(rxode2::etTrans(.d, .m))
      .want <- any(.et$TIME == 1 & .et$EVID == 0)
      expect_identical(.nmObsRows(.d)[2], .want, label = paste(names(.r), .r, collapse = " "))
    }
  })

  test_that(".nmAssertIds() and .nmDropNoObsSubjects() match the FOCEi data step", {
    .d <- data.frame(ID = c(1, 1, 2, 2, 3), TIME = c(0, 1, 0, 1, 0), EVID = c(1, 0, 1, 0, 1), DV = c(NA, 2, NA, NA, NA))
    expect_warning(.k <- .nmDropNoObsSubjects(.d), "IDs without observations dropped: 2 3", fixed = TRUE)
    expect_identical(.k, .d[1:2, ])
    # all subjects have observations: unchanged; none has: kept whole
    expect_silent(expect_identical(.nmDropNoObsSubjects(.d[1:2, ]), .d[1:2, ]))
    expect_silent(expect_identical(.nmDropNoObsSubjects(.d[3:5, ]), .d[3:5, ]))
    .d$ID[2] <- NA
    expect_error(.nmDropNoObsSubjects(.d), "ID is missing on 1 record(s)", fixed = TRUE)
    expect_error(.nmSubjectObs(.d, 1), "ID is missing on 1 record(s)", fixed = TRUE)
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
