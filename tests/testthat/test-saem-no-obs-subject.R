nmTest({
  # Issue #687: SAEM errored ("No data with ID: <id>") for a subject that has a
  # dose but no usable measurement -- all of its observation (EVID==0) rows have
  # DV=NA, which rxode2::etTrans converts to EVID==2 rows, so the subject was not
  # dropped and reached the SAEM kernel with zero observations.  FOCEI tolerated
  # the same subject.  The shared preprocessor (.foceiPreProcessData) now drops
  # subjects that have no usable observation -- emitting the same
  # "IDs without observations dropped" warning rxode2 uses for dose-only
  # subjects -- and the shared table builder (addTable) re-inserts their rows
  # with a population PRED (solved at eta=0) and NA individual columns, so every
  # estimation method behaves like FOCEI.
  .noObsData <- function() {
    d <- nlmixr2data::theo_sd[nlmixr2data::theo_sd$ID %in% 1:3, ]
    # subject 3: dose present, but no measurement (all observations NA)
    d$DV[d$ID == 3 & d$EVID == 0] <- NA_real_
    d
  }

  test_that("SAEM tolerates a dosing subject with no observations like FOCEI (#687)", {
    d <- .noObsData()

    fitSaem <- suppressWarnings(suppressMessages(
      .nlmixr(one.compartment, d, est = "saem", control = saemControlFast)))

    # 1. SAEM no longer errors; it returns a fit
    expect_true(inherits(fitSaem, "nlmixr2FitData"))

    # 2. the dropped subject is reported in the run info (model printing)
    expect_equal(
      sum(grepl("IDs without observations dropped: 3", fitSaem$runInfo, fixed = TRUE)),
      1L)

    # 3. subject 3 is present in the output: a population PRED is filled, but the
    #    measurement and every individual column are NA
    .dfS <- as.data.frame(fitSaem)
    expect_equal(nrow(.dfS), 33L)             # 2 estimated x 11 + 11 re-inserted
    .s3 <- .dfS[as.integer(as.character(.dfS$ID)) == 3, ]
    expect_equal(nrow(.s3), 11L)
    expect_true(all(is.na(.s3$DV)))
    expect_true(all(!is.na(.s3$PRED)))        # population PRED computed
    for (.col in c("IPRED", "eta.ka", "eta.cl", "eta.v",
                   "RES", "IRES", "IWRES")) {
      if (.col %in% names(.s3)) {
        expect_true(all(is.na(.s3[[.col]])), info = .col)
      }
    }

    # 4. estimation used only the two subjects that have data
    expect_setequal(as.integer(as.character(fitSaem$eta$ID)), c(1L, 2L))
  })

  test_that("FOCEI matches SAEM for a dosing subject with no observations (#687)", {
    d <- .noObsData()
    fitFocei <- suppressWarnings(suppressMessages(
      .nlmixr(one.compartment, d, est = "focei", control = foceiControlFast)))
    expect_true(inherits(fitFocei, "nlmixr2FitData"))
    expect_equal(
      sum(grepl("IDs without observations dropped: 3", fitFocei$runInfo, fixed = TRUE)),
      1L)
    .dfF <- as.data.frame(fitFocei)
    .f3 <- .dfF[as.integer(as.character(.dfF$ID)) == 3, ]
    expect_equal(nrow(.f3), 11L)
    expect_true(all(is.na(.f3$DV)))
    expect_true(all(!is.na(.f3$PRED)))        # population PRED computed
    expect_true(all(is.na(.f3$IPRED)))        # individual columns NA
  })

  test_that("multiple no-observation subjects are each dropped and re-inserted (#687)", {
    d <- nlmixr2data::theo_sd[nlmixr2data::theo_sd$ID %in% 1:4, ]
    # subjects 3 and 4: a dose but no measurement
    d$DV[d$ID %in% c(3, 4) & d$EVID == 0] <- NA_real_
    fit <- suppressWarnings(suppressMessages(
      .nlmixr(one.compartment, d, est = "saem", control = saemControlFast)))
    expect_true(inherits(fit, "nlmixr2FitData"))
    # both dropped subjects reported in one message
    expect_equal(
      sum(grepl("IDs without observations dropped: 3 4", fit$runInfo, fixed = TRUE)),
      1L)
    .df <- as.data.frame(fit)
    expect_setequal(as.integer(as.character(unique(.df$ID))), 1:4)
    expect_setequal(as.integer(as.character(fit$eta$ID)), c(1L, 2L))
    for (.id in c(3L, 4L)) {
      .s <- .df[as.integer(as.character(.df$ID)) == .id, ]
      expect_equal(nrow(.s), 11L)
      expect_true(all(is.na(.s$DV)), info = .id)
      expect_true(all(!is.na(.s$PRED)), info = .id)   # population PRED computed
      expect_true(all(is.na(.s$IPRED)), info = .id)   # individual columns NA
    }
  })

  test_that("a fit with all subjects observed is unaffected by the no-obs drop (#687)", {
    d <- nlmixr2data::theo_sd[nlmixr2data::theo_sd$ID %in% 1:3, ]
    fit <- suppressWarnings(suppressMessages(
      .nlmixr(one.compartment, d, est = "saem", control = saemControlFast)))
    expect_true(inherits(fit, "nlmixr2FitData"))
    # no subject dropped -> no such warning, all three subjects present
    expect_equal(
      sum(grepl("IDs without observations dropped", fit$runInfo, fixed = TRUE)),
      0L)
    expect_setequal(as.integer(as.character(unique(fit$ID))), 1:3)
  })
})
