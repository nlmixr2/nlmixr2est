# The covariance-warning sites (functions that call .covWarn() or .covClip()),
# driven with every label a fit really carries: a label is never what gets
# clipped, so each warning is one whole line of at most 74 columns.  Clipping is
# for input no label produces (a long condition message), which the sites with
# such input are driven with too.  The sites in other files are driven there.

test_that("the real labels are the ones a fit carries", {
  .labels <- .covRealLabels()
  expect_true(all(
    c("r,s", "r", "s", "analytic", "linFim", "|r|,s+", "r+,|s|", "|r|,|s|", "|linFim|", "r,s (full)", "|r|,|s| (full)") %in%
      .labels
  ))
  expect_false(anyDuplicated(.labels) > 0L)
  expect_false(any(c("", NA_character_) %in% .labels))
  expect_true(all(.covRealLabels(decorated = FALSE, full = FALSE) %in% .labels))
  expect_false(any(vapply(.covRealLabels(full = FALSE), .covIsFull, logical(1))))
  expect_false(any(grepl("|", .covRealLabels(decorated = FALSE), fixed = TRUE)))
})

test_that(".covInstall() says what it installed instead of the request, whole, for every label", {
  .labels <- .covRealLabels()
  .got <- character()
  .want <- character()
  .covSite(".covInstall", {
    for (.have in .labels) {
      for (.asked in .labels) {
        if (.covSameName(.covBaseName(.asked), .covBaseName(.have))) {
          next
        }
        .got <- c(.got, capture_warnings(.covInstall(.covSiteEnv(), .covSitePd(), .have, what = .asked)))
        .want <- c(.want, sprintf("\"%s\" installed instead of the requested \"%s\"", .have, .asked))
      }
    }
  })
  expect_gt(length(.want), 1000L)
  expect_identical(.got, .want)
  expect_true(.covSiteWhole(.got))
})

test_that(".covRejectWarn() is one whole line for every label and reason", {
  .labels <- .covRealLabels()
  .reasons <- c(
    "is not finite",
    "is not symmetric",
    "is not positive definite",
    "could not be computed",
    "needs a positive-definite R",
    "needs a positive-definite S"
  )
  .got <- character()
  .want <- character()
  .covSite(".covRejectWarn", {
    for (.what in .labels) {
      for (.kept in c(NA_character_, .labels)) {
        for (.why in .reasons) {
          .e <- new.env(parent = emptyenv())
          if (!is.na(.kept)) {
            .e$cov <- .covSitePd()
            .e$covMethod <- .kept
          }
          .got <- c(.got, capture_warnings(.covRejectWarn(.e, .what, .why)))
          .want <- c(
            .want,
            sprintf(
              "\"%s\" %s; %s",
              .what,
              .why,
              if (is.na(.kept)) "none installed" else sprintf("kept \"%s\"", .kept)
            )
          )
        }
      }
    }
  })
  expect_identical(.got, .want)
  expect_true(.covSiteWhole(.got))
})

test_that(".covRepairWarn() is one whole line for every repaired label", {
  # a repaired covariance that is not a finite-difference one (SAEM's "|linFim|",
  # "|Ha|") is installed without the " (full)" scope suffix
  .repaired <- Filter(
    function(l) !identical(.covUnrepaired(l), l) && !nzchar(.covFdType(l)),
    .covRealLabels(full = FALSE)
  )
  expect_gt(length(.repaired), 1L)
  .got <- character()
  .covSite(".covRepairWarn", {
    for (.l in .repaired) {
      .got <- c(.got, capture_warnings(.covRepairWarn(.covUnrepaired(.l), .l)))
    }
  })
  expect_identical(
    .got,
    sprintf("\"%s\" is not positive definite; installed as \"%s\" (sqrtm)", .covUnrepaired(.repaired), .repaired)
  )
  expect_true(.covSiteWhole(.got))
})

test_that(".covNotUsableMsg() is one whole line for every request, reason and installed label", {
  .asked <- .covRealLabels(decorated = FALSE, full = FALSE)
  .tails <- c("no covariance", sprintf("installed \"%s\"", .covRealLabels(full = FALSE)))
  .reasons <- c("R not PD", "S not PD", "R not PD; S not PD", "S failed", "sandwich doubtful", "all variances < 1e-7")
  .grid <- expand.grid(request = .asked, why = .reasons, tail = .tails, stringsAsFactors = FALSE)
  .got <- .covSite(
    ".covNotUsableMsg",
    mapply(.covNotUsableMsg, .grid$request, .grid$why, .grid$tail, USE.NAMES = FALSE)
  )
  expect_identical(.got, sprintf("\"%s\" not usable (%s); %s", .grid$request, .grid$why, .grid$tail))
  expect_true(.covSiteWhole(.got))
})

test_that(".foceiAnalyticClip() keeps a short reason and clips a long one inside its width", {
  .covSite(".foceiAnalyticClip", {
    expect_identical(.foceiAnalyticClip("an error"), "an error")
    expect_identical(.foceiAnalyticClip(strrep("x", 80)), paste0(strrep("x", 51), "..."))
  })
})

test_that("the analytic step error is one line whatever the message", {
  withr::defer(.foceiAnalyticWhyReset())
  .msgs <- c("no solve", paste0("rxSolve: ", strrep("m", 80)), "first line\n   second line")
  .got <- character()
  .covSite(".foceiAnalyticErrWarn", {
    for (.site in 1:2) {
      for (.m in .msgs) {
        .got <- c(.got, capture_warnings(.foceiAnalyticErrWarn(.site)(simpleError(.m))))
      }
    }
  })
  expect_identical(
    .got,
    c(
      "analytic err1: no solve",
      paste0("analytic err1: rxSolve: ", strrep("m", 42), "..."),
      "analytic err1: first line second line",
      "analytic err2: no solve",
      paste0("analytic err2: rxSolve: ", strrep("m", 42), "..."),
      "analytic err2: first line second line"
    )
  )
  expect_true(all(cli::ansi_nchar(.got, type = "width") <= 74L))
})

test_that(".foceiWarnConditionalCov() names every finite-difference label whole", {
  .fd <- Filter(function(l) nzchar(.covFdType(l)), .covRealLabels())
  expect_gt(length(.fd), 10L)
  .got <- character()
  .covSite(".foceiWarnConditionalCov", {
    for (.l in .fd) {
      .e <- new.env(parent = emptyenv())
      .e$etaObf <- 1
      .e$cov <- diag(2)
      .e$covMethod <- .l
      .got <- c(.got, capture_warnings(.foceiWarnConditionalCov(.e, list(maxInnerIterations = 0L))))
    }
  })
  expect_length(.got, length(.fd))
  expect_identical(startsWith(.got, sprintf("\"%s\" ", .fd)), rep(TRUE, length(.fd)))
  expect_true(.covSiteWhole(.got))
})

test_that(".setCovFd() names the label it installed and the one asked for, whole", {
  .fd <- Filter(function(l) nzchar(.covFdType(l)), .covRealLabels())
  .cur <- NULL
  local_mocked_bindings(
    .setCovRefit = function(...) list(covMethod = .cur, cov = .covSitePd(), parFixedDf = NULL, parFixed = NULL),
    .covInstall = function(...) TRUE,
    .covCacheDrop = function(...) NULL
  )
  .got <- character()
  .want <- character()
  .covSite(".setCovFd", {
    for (.asked in .fd) {
      for (.have in .fd) {
        if (.covSameName(.asked, .have) || !identical(.covIsFull(.asked), .covIsFull(.have))) {
          next
        }
        .cur <- .have
        .ctl <- list(covFallback = stats::setNames(list(c("r", "s", "r,s")), .covBaseName(.asked)))
        .got <- c(.got, capture_warnings(.setCovFd(NULL, new.env(), .asked, .covBaseName(.asked), control = .ctl)))
        .want <- c(.want, sprintf("\"%s\" installed instead of the requested \"%s\"", .have, .asked))
      }
    }
  })
  expect_gt(length(.want), 20L)
  expect_identical(.got, .want)
  expect_true(.covSiteWhole(.got))
})

test_that(".foceiDropMalformedCov() drops a malformed covariance with one whole line", {
  .nm <- c("tka", "tcl")
  .bad <- matrix(c(1, NA, NA, 1), 2, dimnames = list(.nm, .nm))
  .e <- new.env(parent = emptyenv())
  .e$cov <- .bad
  .w <- .covSite(".foceiDropMalformedCov", capture_warnings(suppressMessages(.foceiDropMalformedCov(.e))))
  expect_identical(.w, "covariance not in proper form, can access value in $covDebug")
  expect_null(.e$cov)
  expect_identical(.e$covDebug, .bad)
  # whole ill-identified rows and columns, and a proper matrix, are kept
  .part <- matrix(c(1, NA, NA, NA), 2, dimnames = list(.nm, .nm))
  for (.ok in list(.part, .covSitePd(.nm))) {
    .e <- new.env(parent = emptyenv())
    .e$cov <- .ok
    expect_silent(.foceiDropMalformedCov(.e))
    expect_identical(.e$cov, .ok)
  }
  # no covariance, nothing to check
  expect_silent(.foceiDropMalformedCov(new.env(parent = emptyenv())))
})
