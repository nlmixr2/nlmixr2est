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
    # as etTrans() reads them: EVID and MDV truncated, NaN not missing, and a
    # negative or "-name" CMT turns a compartment off
    .d <- data.frame(
      TIME = c(1, NaN, 1, 1, 1, 1),
      EVID = c(0.7, 0, 1.2, 0, 0, 0),
      DV = c(4, 4, 4, NaN, 4, 4),
      MDV = c(0, 0, 0, 0, 1.5, 0),
      CMT = c(2, 2, 2, 2, 2, -2)
    )
    expect_identical(.nmObsRows(.d), c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE))
    expect_identical(.nmObsRows(data.frame(TIME = 1, EVID = 0, DV = 4, CMT = c("center", "-center"))), c(TRUE, FALSE))
    # the columns are matched exactly: DVID is not DV
    expect_error(.nmObsRows(data.frame(TIME = 0, EVID = 0, DVID = 1)), "needs DV, EVID and TIME columns")
    expect_error(.nmObsRows(data.frame(TIME = 0, DV = 1)), "needs DV, EVID and TIME columns")
    expect_error(.nmObsRows(data.frame(EVID = 0, DV = 1)), "needs DV, EVID and TIME columns")
  })

  # The message .foceiPreProcessData() stops with on `d`, or NULL; `env` keeps
  # what it builds
  .utFoceiError <- function(d, ui, env = new.env(parent = emptyenv())) {
    tryCatch(
      {
        suppressWarnings(.foceiPreProcessData(d, env, ui))
        NULL
      },
      error = function(e) conditionMessage(e)
    )
  }

  # FOCEi's verdict on the second of three records (a dose, the record, and an
  # observation at t = 2): TRUE when its data ($dataSav, after its own defaults
  # and etTrans()) has two EVID 0 records.  When etTrans() rejects the data, its
  # message, which must be one of the reasons no data can avoid: a dose with a
  # missing AMT, EVID 1 or 4 with an MDV of 0 or missing, or (rxode2 after 5.1.8)
  # a missing EVID on a dose record
  .utFoceiReasons <- paste0(
    "'amt' value NA for dose event|'mdv' cannot be 0 when 'evid'=[14]",
    if (utils::packageVersion("rxode2") > "5.1.8") {
      "|'evid' is missing on a record with dose information"
    }
  )
  .utFoceiObs <- function(d, ui) {
    .e <- new.env(parent = emptyenv())
    .err <- .utFoceiError(d, ui, .e)
    if (!is.null(.err)) {
      expect_match(.err, .utFoceiReasons, label = paste(names(d), d[2, ], collapse = " "))
      return(.err)
    }
    sum(.e$dataSav$EVID == 0) == 2L
  }

  # a dose, the record under test and an observation at t = 2; a NULL evid, amt
  # or mdv leaves that column out
  .utObsRecord <- function(time, dv, cmt, evid, amt, mdv) {
    .d <- data.frame(ID = 1L, TIME = c(0, time, 2), DV = c(NA, dv, 3), CMT = c(1, cmt, 2))
    if (!is.null(evid)) {
      .d$EVID <- c(1, evid, 0)
    }
    if (!is.null(amt)) {
      .d$AMT <- c(100, amt, 0)
    }
    if (!is.null(mdv)) {
      .d$MDV <- c(1, mdv, 0)
    }
    .d
  }

  test_that(".nmObsRows() and .nmWithEvid() agree with FOCEi's translated data", {
    skip_on_cran()
    .ui <- rxode2::assertRxUi(function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.cl ~ 0.3
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    })
    # -1 (EVID, MDV) and Inf (AMT): no such column
    .g <- expand.grid(
      evid = c(0, 1, 2, 3, 4, 101, NA, -1),
      dv = c(4, NA, NaN),
      mdv = c(-1, 0, 1, 2, NA),
      amt = c(0, NA, 50, -50, Inf),
      time = c(1, NA),
      cmt = c(2, -2)
    )
    .n <- 0L
    .nNew <- 0L
    for (.k in seq_len(nrow(.g))) {
      .r <- .g[.k, ]
      .d <- .utObsRecord(
        .r$time,
        .r$dv,
        .r$cmt,
        if (.r$evid %in% -1) NULL else .r$evid,
        if (.r$amt %in% Inf) NULL else .r$amt,
        if (.r$mdv %in% -1) NULL else .r$mdv
      )
      .want <- .utFoceiObs(.d, .ui)
      if (is.character(.want)) {
        .nNew <- .nNew + grepl("'evid' is missing on a record", .want, fixed = TRUE)
        next
      }
      .n <- .n + 1L
      expect_identical(.nmObsRows(.nmWithEvid(.d))[2], .want, label = paste(names(.r), .r, collapse = " "))
    }
    # 2400 combinations, of which etTrans() rejects 174 for the first two reasons;
    # a later rxode2 also rejects the missing-EVID ones with a nonzero AMT
    expect_identical(.n + .nNew, 2226L)
    # values the grid leaves out: NaN TIME, AMT, EVID and MDV, fractional EVID and
    # MDV, and a "-name" CMT
    .base <- list(time = 1, dv = 4, cmt = 2, evid = 0, amt = 0, mdv = 0)
    .special <- list(
      list(time = NaN),
      list(dv = NA, amt = NaN),
      list(evid = NaN),
      list(mdv = NaN),
      list(evid = 0.7),
      list(evid = 1.2, mdv = 1),
      list(mdv = 0.5),
      list(mdv = 1.5),
      list(evid = NULL, mdv = 0.5, amt = 50),
      list(cmt = "-center")
    )
    for (.v in .special) {
      .a <- utils::modifyList(.base, .v, keep.null = TRUE)
      .d <- do.call(.utObsRecord, .a)
      if (identical(.a$cmt, "-center")) {
        .d$CMT <- c("depot", "-center", "center")
      }
      expect_identical(
        .nmObsRows(.nmWithEvid(.d))[2],
        .utFoceiObs(.d, .ui),
        label = paste(names(.v), .v, collapse = " ")
      )
    }
  })

  test_that(".nmAssertEventColumns() rejects the column types FOCEi rejects, with its message", {
    skip_on_cran()
    .ui <- rxode2::assertRxUi(function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.cl ~ 0.3
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    })
    .base <- data.frame(
      ID = 1L,
      TIME = c(0, 1, 2),
      DV = c(NA, 4, 3),
      AMT = c(100, 0, 0),
      EVID = c(1, 0, 0),
      CMT = c(1, 2, 2),
      MDV = c(1, 0, 0)
    )
    .types <- list(
      character = as.character,
      factor = factor,
      logical = as.logical,
      naOnly = function(x) rep(NA, length(x)),
      Date = function(x) as.Date(x, origin = "1970-01-01"),
      POSIXct = function(x) as.POSIXct(x, origin = "1970-01-01", tz = "UTC")
    )
    # with and without the EVID and AMT columns, which change the MDV message
    .layouts <- list(all = character(0), noEvid = "EVID", noEvidAmt = c("EVID", "AMT"))
    .known <- "needs to be an integer|needs to be a number|infinite times are not allowed"
    .n <- 0L
    for (.lay in names(.layouts)) {
      for (.col in c("EVID", "MDV", "AMT", "TIME", "DV")) {
        if (.col %in% .layouts[[.lay]]) {
          next
        }
        for (.ty in names(.types)) {
          .d <- .base[!(names(.base) %in% .layouts[[.lay]])]
          .d[[.col]] <- .types[[.ty]](.base[[.col]])
          .lab <- paste(.lay, .col, .ty)
          .f <- .utFoceiError(.d, .ui)
          .h <- tryCatch(
            {
              .nmAssertEventColumns(.d)
              NULL
            },
            error = function(e) conditionMessage(e)
          )
          if (is.null(.h)) {
            # FOCEi may stop for another reason, but not for one this check covers
            expect_false(!is.null(.f) && grepl(.known, .f), label = paste(.lab, .f))
          } else {
            expect_identical(.h, .f, label = .lab)
            .n <- .n + 1L
          }
        }
      }
    }
    # the rejections are the character EVID, MDV and AMT, in each layout that has
    # them (3 + 2 + 1); a factor, a logical, a Date or a POSIXct is
    # stored as a number, as etTrans() checks it
    expect_identical(.n, 6L)
    # infinite times: the first in the data as given, in several subjects out of
    # order, and with no ID column
    .inf <- data.frame(ID = c(2L, 2L, 1L, 1L), TIME = c(1, 0, -Inf, Inf), DV = c(1, 2, 3, 4), EVID = 0)
    expect_identical(
      tryCatch(.nmAssertEventColumns(.inf), error = conditionMessage),
      .utFoceiError(.inf, .ui)
    )
    .inf$ID <- NULL
    expect_identical(
      tryCatch(.nmAssertEventColumns(.inf), error = conditionMessage),
      .utFoceiError(.inf, .ui)
    )
    # a missing ID beside each of those: FOCEi checks any ID but a double first,
    # and a double one record by record with TIME, after the event columns
    .ids <- list(double = as.double, integer = as.integer, character = as.character, factor = factor)
    .m <- 0L
    for (.it in names(.ids)) {
      .cases <- list()
      for (.lay in c("all", "noEvid")) {
        for (.col in c("EVID", "MDV", "AMT")) {
          if (.col %in% .layouts[[.lay]]) {
            next
          }
          .d <- .base[!(names(.base) %in% .layouts[[.lay]])]
          .d[[.col]] <- as.character(.d[[.col]])
          .d$ID <- .ids[[.it]](c(NA, 1, 1))
          .cases[[paste(.lay, .col)]] <- .d
        }
      }
      .d <- .base
      .d$ID <- .ids[[.it]](c(NA, 1, 1))
      .d$TIME[3] <- Inf
      .cases[["ID before TIME"]] <- .d
      .d <- .base
      .d$ID <- .ids[[.it]](c(1, 1, NA))
      .d$TIME[1] <- -Inf
      .cases[["TIME before ID"]] <- .d
      for (.nm in names(.cases)) {
        .f <- .utFoceiError(.cases[[.nm]], .ui)
        expect_false(is.null(.f), label = paste(.it, .nm))
        expect_identical(
          tryCatch(.nmAssertEventColumns(.cases[[.nm]]), error = conditionMessage),
          .f,
          label = paste(.it, .nm)
        )
        .m <- .m + 1L
      }
    }
    expect_identical(.m, 28L)
    # a logical ID: etTrans() cannot read one at all (a type-conversion error),
    # with or without a missing value; the check stops first and names it
    for (.l in list(c(TRUE, TRUE, TRUE), c(NA, TRUE, TRUE))) {
      .d <- .base
      .d$ID <- .l
      expect_false(is.null(.utFoceiError(.d, .ui)))
      expect_identical(tryCatch(.nmAssertEventColumns(.d), error = conditionMessage), "'id' cannot be logical")
    }
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
    # a missing ID, with etTrans()'s message for a double ID and for any other
    .d$ID[2] <- NA
    expect_error(.nmDropNoObsSubjects(.d), "^'id' cannot be 'NA'$")
    expect_error(.nmSubjectObs(.d, 1), "^'id' cannot be 'NA'$")
    .d$ID <- as.integer(.d$ID)
    expect_error(.nmDropNoObsSubjects(.d), "^'id' cannot have NA values$")
    expect_error(.nmSubjectObs(.d, 1), "^'id' cannot have NA values$")
    # with no TIME column the record check still finds a missing double ID, and a
    # missing ID of any other type is checked first anyway
    expect_error(.nmAssertEventColumns(data.frame(ID = c(NA, 1), DV = 1, EVID = 0)), "^'id' cannot be 'NA'$")
    expect_error(.nmAssertEventColumns(data.frame(ID = c(NA, 1L), DV = 1, EVID = 0)), "^'id' cannot have NA values$")
    expect_error(.nmWithEvid(data.frame(ID = c(NA, 1), DV = 1)), "^'id' cannot be 'NA'$")
    expect_error(.nmWithEvid(data.frame(ID = c(NA, 1L), DV = 1)), "^'id' cannot have NA values$")
    # a missing double ID beside a character EVID: the EVID first, as FOCEi checks
    # a double ID with the records, after the event columns
    .e <- .d
    .e$ID <- as.double(.e$ID)
    .e$EVID <- as.character(.e$EVID)
    expect_error(.nmDropNoObsSubjects(.e), "^event id \\('evid'\\) needs to be an integer$")
    expect_error(.nmSubjectObs(.e, 1), "^event id \\('evid'\\) needs to be an integer$")
    # the columns are matched exactly: IDENT is not ID
    .i <- data.frame(IDENT = 1, TIME = 0, EVID = 0, DV = 1)
    expect_error(.nmDropNoObsSubjects(.i), "the data needs an ID column", fixed = TRUE)
    expect_error(.nmSubjectObs(.i, 1), "the data needs an ID column", fixed = TRUE)
    expect_error(.nmSubjectObs(data.frame(ID = 1, EVID = 0, DV = 1), 1), "the data needs a TIME column", fixed = TRUE)
  })

  test_that(".nmObsRows() and .nmWithEvid() stop on the columns FOCEi rejects, with its message", {
    .d <- data.frame(ID = 1, TIME = c(0, 1), EVID = c(1, 0), DV = c(NA, 1), AMT = c(100, 0), MDV = c(1, 0))
    .e <- .d
    .e$EVID <- as.character(.e$EVID)
    expect_error(.nmObsRows(.e), "event id ('evid') needs to be an integer", fixed = TRUE)
    .m <- .d
    .m$MDV <- as.character(.m$MDV)
    expect_error(.nmObsRows(.m), "missing dependent variable ('mdv') needs to be an integer", fixed = TRUE)
    # with no EVID column but an AMT, the MDV column stands for the EVID
    expect_error(.nmWithEvid(.m[names(.m) != "EVID"]), "^missing DV \\('mdv'\\) needs to be an integer$")
    .n <- .m[!(names(.m) %in% c("EVID", "AMT"))]
    expect_error(.nmWithEvid(.n), "^missing dependent variable \\('mdv'\\) needs to be an integer$")
    .a <- .d
    .a$AMT <- as.character(.a$AMT)
    expect_error(.nmObsRows(.a), "amount ('amt') needs to be a number", fixed = TRUE)
    expect_error(.nmWithEvid(.a[names(.a) != "EVID"]), "amount ('amt') needs to be a number", fixed = TRUE)
    for (.inf in c(Inf, -Inf)) {
      .t <- .d
      .t$TIME[2] <- .inf
      expect_error(.nmObsRows(.t), "infinite times are not allowed (id: 1, row: 2)", fixed = TRUE)
      expect_error(.nmWithEvid(.t), "infinite times are not allowed (id: 1, row: 2)", fixed = TRUE)
      # without an ID column the subject is 1, as FOCEi makes it
      expect_error(.nmObsRows(.t[names(.t) != "ID"]), "infinite times are not allowed (id: 1, row: 2)", fixed = TRUE)
    }
    # the first infinite time in the data as given
    .o <- data.frame(ID = c(2, 2, 1, 1), TIME = c(1, 0, -Inf, Inf), EVID = 0, DV = 1)
    expect_error(.nmObsRows(.o), "infinite times are not allowed (id: 1, row: 3)", fixed = TRUE)
    # a missing double ID ahead of the first infinite time is reported, as FOCEi
    # checks the two record by record
    .o$ID[1] <- NA
    expect_error(.nmObsRows(.o), "^'id' cannot be 'NA'$")
    # what FOCEi reads as numbers: a character TIME or DV, a logical or factor EVID
    .c <- .d
    .c$TIME <- c("0", "1")
    .c$DV <- c(".", "1")
    expect_identical(.nmObsRows(.c), c(FALSE, TRUE))
    .c$EVID <- c(TRUE, FALSE)
    expect_identical(.nmObsRows(.c), c(FALSE, TRUE))
  })

  test_that(".nmWithEvid() gives data without an EVID column FOCEi's", {
    # no AMT column either: every record EVID 0, and an MDV column stays an MDV
    expect_identical(.nmWithEvid(data.frame(DV = 1:2))$EVID, c(0L, 0L))
    .m <- .nmWithEvid(data.frame(DV = 1:2, MDV = c(1, 0)))
    expect_identical(.m$EVID, c(0L, 0L))
    expect_identical(.m$MDV, c(1, 0))
    # a nonzero AMT (NaN too, as etTrans() reads it) is a dose
    expect_identical(.nmWithEvid(data.frame(AMT = c(0, NA, 100, -1, NaN)))$EVID, c(0L, 0L, 1L, 1L, 1L))
    # with AMT, an MDV column is the EVID and is dropped: a missing MDV is a dose,
    # MDV 1 with a zero AMT is EVID 2, and other values are truncated
    .m <- .nmWithEvid(data.frame(AMT = c(100, 0, 50, 0, 0, 0), MDV = c(1, 1, 0, NA, 2, 0.5)))
    expect_identical(.m$EVID, c(1L, 2L, 0L, 1L, 2L, 0L))
    expect_null(.m$MDV)
    # data with an EVID column is unchanged
    .k <- data.frame(EVID = c(0, 1), AMT = c(0, 1), MDV = c(0, 1))
    expect_identical(.nmWithEvid(.k), .k)
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
    # a reset starts a new episode at an earlier clock time; episodes keep their order
    .r <- data.frame(
      ID = 1,
      TIME = c(0, 1, 0, 1),
      EVID = c(1, 0, 4, 0),
      AMT = c(100, 0, 200, 0),
      DV = c(NA, 3, NA, 6)
    )
    .s <- .nmSubjectObs(.r, 1)[[1]]
    expect_identical(.s$ev, .r)
    expect_identical(.s$times, c(1, 1))
    expect_identical(.s$y, c(3, 6))
    .s <- .nmSubjectObs(.r[c(2, 1, 3, 4), ], 1)[[1]]
    expect_identical(.s$ev, .r)
    expect_identical(.s$y, c(3, 6))
  })
})
