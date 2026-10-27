# Solving a finished fit with inter-occasion variability again: rxSolve(),
# predict() at both levels and augPred().  The individual solves read each
# occasion eta from the fit's estimates, so at the fit's records they give
# its IPRED; the population solves set every eta to zero and give its PRED.

# The fits shared by the tests below, each made once per run
.iovPredictCache <- new.env(parent = emptyenv())

# The cached `.iovOdeFit(...)` named `key`
.iovPredictFit <- function(key, ...) {
  if (!exists(key, envir = .iovPredictCache, inherits = FALSE)) {
    assign(key, .iovOdeFit(...), envir = .iovPredictCache)
  }
  get(key, envir = .iovPredictCache, inherits = FALSE)
}

# The fits whose augPred() is checked: estimation method and how the occasion
# changes
.iovPredictAugCases <- list(
  list(key = "focei", est = "focei"),
  list(key = "saem", est = "saem"),
  list(key = "focei-subject", est = "focei", occ = "subject"),
  list(key = "focei-nocb", est = "focei", interpolation = "nocb"),
  list(key = "focei-tie-obs-dose", est = "focei", data = .iovOdeTieData("obs-dose")),
  list(key = "focei-tie-dose-obs", est = "focei", data = .iovOdeTieData("dose-obs")),
  list(key = "focei-two-occ", est = "focei", model = .iovOdeTwoOccMod),
  list(key = "focei-id-shifted", est = "focei", data = .iovOdeIdData("shifted")),
  list(key = "focei-id-character", est = "focei", data = .iovOdeIdData("character")),
  list(key = "focei-id-descending", est = "focei", data = .iovOdeIdData("descending")),
  list(key = "focei-id-dose-only", est = "focei", data = .iovOdeIdData("dose-only")),
  list(key = "focei-addl", est = "focei", data = .iovOdeAddlData())
)

# The cached fit of one of .iovPredictAugCases
.iovPredictCaseFit <- function(case) {
  do.call(.iovPredictFit, c(list(case$key), case[setdiff(names(case), "key")]))
}

nmTest({
  test_that("a fit with IOV keeps its nested omega, so it can be solved", {
    for (.est in c("focei", "saem")) {
      .fit <- .iovPredictFit(.est, .est)
      expect_s3_class(.fit$omega, "lotri")
      expect_identical(names(.fit$omega), c("id", "occ"), label = .est)
      .s <- suppressWarnings(suppressMessages(rxode2::rxSolve(.fit)))
      expect_identical(nrow(as.data.frame(.s)), nrow(.fit), label = .est)
    }
  })

  test_that("individual predict() of a fit with IOV uses its occasion etas", {
    for (.est in c("focei", "saem")) {
      .fit <- .iovPredictFit(.est, .est)
      .obs <- as.data.frame(.fit)
      set.seed(1)
      .p1 <- suppressMessages(predict(.fit, level = "individual"))
      set.seed(2)
      .p2 <- suppressMessages(predict(.fit, level = "individual"))
      .pn <- suppressMessages(predict(.fit, newdata = .iovOdeData(), level = "individual"))
      # no eta is simulated, so the seed changes nothing
      expect_identical(.p1$sim, .p2$sim, label = .est)
      expect_identical(.p1$iov.cl, .p2$iov.cl, label = .est)
      expect_equal(.p1$sim, .obs$IPRED, tolerance = 1e-8, label = .est)
      expect_equal(.pn$sim, .obs$IPRED, tolerance = 1e-8, label = .est)
      .iov <- .fit$env$iov$occ
      .k <- match(paste(.obs$ID, ifelse(.obs$TIME < 6, 1L, 2L)), paste(.iov$ID, .iov$occ))
      expect_identical(.p1$iov.cl, .iov$iov.cl[.k], label = .est)
    }
  })

  test_that("individual predict() joins every eta on the subject, in any row order", {
    # newdata in descending subject order: each row still takes its own
    # subject's etas, not those of the subject in its position
    for (.est in c("focei", "saem")) {
      .fit <- .iovPredictFit(.est, .est)
      .obs <- as.data.frame(.fit)
      .d <- .iovOdeIdData("descending")
      .p <- suppressMessages(predict(.fit, newdata = .d, level = "individual"))
      .k <- match(paste(.obs$ID, .obs$TIME), paste(.p$id, .p$time))
      expect_false(anyNA(.k), label = .est)
      expect_equal(.p$sim[.k], .obs$IPRED, tolerance = 1e-8, label = .est)
      expect_identical(
        .p$eta.v[.k],
        .fit$eta$eta.v[match(as.character(.obs$ID), as.character(.fit$eta$ID))],
        label = .est
      )
    }
  })

  test_that("individual predict() of a fit whose occasion starts at an addl dose gives its IPRED", {
    # rxode2 adds the 12 h dose with no occasion (addlKeepsCov = FALSE), and
    # fills the eta column there as the fit's own solve fills the occasion
    .fit <- .iovPredictCaseFit(.iovPredictAugCases[[length(.iovPredictAugCases)]])
    .p <- suppressMessages(predict(.fit, level = "individual"))
    expect_equal(.p$sim, as.data.frame(.fit)$IPRED, tolerance = 1e-8)
  })

  test_that("individual predict() gives an occasion the fit did not see eta 0, with a warning", {
    .fit <- .iovPredictFit("focei", "focei")
    .d <- .iovOdeData()
    .d$occ[.d$ID == 1 & .d$TIME > 12] <- 3L
    .w <- capture_warnings(.p <- suppressMessages(predict(.fit, newdata = .d, level = "individual")))
    expect_identical(.w, "1 subject-occasion pair(s) not in the fit get occasion etas of 0")
    .new <- .p$id == "1" & .p$time > 12
    expect_true(any(.new))
    expect_true(all(.p$iov.cl[.new] == 0))
    expect_false(any(.p$iov.cl[!.new] == 0))
  })

  test_that("population predict() of a fit with IOV gives its PRED", {
    for (.est in c("focei", "saem")) {
      .fit <- .iovPredictFit(.est, .est)
      .p <- suppressMessages(predict(.fit))
      expect_equal(.p$pred, as.data.frame(.fit)$PRED, tolerance = 1e-8, label = .est)
    }
  })

  test_that("augPred() on a fit with IOV gives the fit's predictions at its records", {
    # the occasion of each added time follows the fit's covariate
    # interpolation and, at tied records, their record order; augPred()'s
    # finer grid solves the ODEs again, to its own tolerance
    for (.case in .iovPredictAugCases) {
      .fit <- .iovPredictCaseFit(.case)
      .obs <- as.data.frame(.fit)
      .a <- as.data.frame(suppressWarnings(suppressMessages(augPred(.fit))))
      expect_identical(levels(.a$ind), c("Individual", "Population", "Observed"), label = .case$key)
      .ind <- .a[.a$ind == "Individual", ]
      .pop <- .a[.a$ind == "Population", ]
      .ki <- match(paste(.obs$ID, .obs$TIME), paste(.ind$id, .ind$time))
      .kp <- match(paste(.obs$ID, .obs$TIME), paste(.pop$id, .pop$time))
      expect_false(anyNA(.ki), label = .case$key)
      expect_false(anyNA(.kp), label = .case$key)
      expect_equal(.ind$values[.ki], .obs$IPRED, tolerance = 1e-5, label = .case$key)
      expect_equal(.pop$values[.kp], .obs$PRED, tolerance = 1e-5, label = .case$key)
    }
  })

  test_that("augPred() gives a time before a subject's first record its first occasion", {
    .fit <- .iovPredictFit("focei", "focei")
    .ev <- .augPredIovEvents(.fit, .augPredExpandData(.fit, minimum = -2))
    .early <- .ev[.ev$TIME < 0, ]
    expect_true(nrow(.early) > 0L)
    .iov <- .fit$env$iov$occ
    .k <- match(paste(as.character(.early$ID), 1L), paste(.iov$ID, .iov$occ))
    expect_identical(.early$iov.cl, .iov$iov.cl[.k])
  })

  test_that("under midpoint interpolation, IOV predictions are the fit's until an occasion changes", {
    # the fit compares the interpolated occasion with each occasion's value,
    # so its occasion etas are 0 between records of two occasions; the
    # solves of the user's model interpolate the etas there instead
    .fit <- .iovPredictFit("focei-midpoint", "focei", interpolation = "midpoint")
    .obs <- as.data.frame(.fit)
    .before <- .obs$TIME < 5
    .p <- suppressMessages(predict(.fit, level = "individual"))
    expect_equal(.p$sim[.before], .obs$IPRED[.before], tolerance = 1e-8)
    .a <- as.data.frame(suppressWarnings(suppressMessages(augPred(.fit))))
    .ind <- .a[.a$ind == "Individual", ]
    .k <- match(paste(.obs$ID, .obs$TIME)[.before], paste(.ind$id, .ind$time))
    expect_equal(.ind$values[.k], .obs$IPRED[.before], tolerance = 1e-5)
    expect_equal(suppressMessages(predict(.fit))$pred, .obs$PRED, tolerance = 1e-8)
  })
})
