nmTest({
  skip_if_not(.nlmixr2estRxode2HasTte(), "rxode2 without tte() endpoints")
  .dir <- system.file("validation", "tte", package = "nlmixr2est")
  skip_if(.dir == "", "tte validation scenarios not installed")
  .env <- new.env()
  suppressMessages(sys.source(file.path(.dir, "scenarios.R"), envir = .env))

  test_that("every tte validation scenario has a data set nlmixr2 can read", {
    .all <- .env$tteValidationScenarios()
    expect_true(length(.all) >= 10)
    for (.name in names(.all)) {
      .s <- .all[[.name]]
      .ep <- suppressMessages(.env$.tteEndpoints(.s))
      .d <- utils::read.csv(file.path(.dir, "data", paste0(.name, ".csv")))
      expect_equal(length(unique(.d$ID)), .s$n, info = .name)
      expect_true(all(c("ID", "TIME", "EVID", "AMT", "CMT", "DV", "MDV", "DVID") %in% names(.d)), info = .name)
      expect_false(anyNA(.d), info = .name)
      expect_equal(.d$MDV, as.integer(.d$EVID != 0), info = .name)
      expect_true(all(diff(.d$TIME)[diff(.d$ID) == 0] >= 0), info = .name)
      .n <- suppressMessages(.env$tteValidationData(.name, .dir))
      expect_true(all(.n$cmt[.n$evid == 0] %in% .ep), info = .name)
      # each tte endpoint starts with a DV = 0 record and scores 0/1 after it
      .tte <- .ep[suppressMessages(rxode2::rxode2(.s$model))$predDf$distribution == "tte"]
      for (.e in .tte) {
        .o <- .n[.n$evid == 0 & .n$cmt == .e, ]
        expect_true(all(.o$dv %in% c(0, 1)), info = paste(.name, .e))
        expect_true(all(.o$dv[!duplicated(.o$id)] == 0), info = paste(.name, .e))
      }
      # the starting values name parameters of the model
      .ini <- suppressMessages(rxode2::rxode2(.s$model))$iniDf$name
      expect_true(all(names(.s$start) %in% .ini), info = .name)
    }
  })

  test_that("the NONMEM columns hold the value in force before each record", {
    .d <- utils::read.csv(file.path(.dir, "data", "tvCov.csv"))
    .prev <- c(NA, .d$X[-nrow(.d)])
    .first <- !duplicated(.d$ID)
    expect_equal(.d$XNM[!.first], .prev[!.first])
    expect_equal(.d$XNM[.first], .d$X[.first])
    expect_equal(.d$X, as.numeric(.d$TIME >= .d$TSW))
    .d <- utils::read.csv(file.path(.dir, "data", "gapTime.csv"))
    for (.id in unique(.d$ID)[1:20]) {
      .s <- .d[.d$ID == .id, ]
      .ev <- .s$TIME[.s$DV == 1]
      .want <- vapply(.s$TIME, function(t) max(c(0, .ev[.ev < t])), numeric(1))
      expect_equal(.s$TLAST, .want)
    }
  })

  test_that("a dropout right-censors the other endpoint", {
    .d <- utils::read.csv(file.path(.dir, "data", "dropout.csv"))
    .drop <- .d[.d$DVID == 2 & .d$DV == 1, ]
    expect_true(nrow(.drop) > 0)
    for (.i in seq_len(nrow(.drop))) {
      .ev <- .d[.d$ID == .drop$ID[.i] & .d$DVID == 1, ]
      # either the event came first, or `ev` is censored at the dropout
      expect_true(any(.ev$DV == 1) || .drop$TIME[.i] %in% .ev$TIME[.ev$DV == 0])
      expect_true(all(.d$TIME[.d$ID == .drop$ID[.i]] <= .drop$TIME[.i]))
    }
  })
})
