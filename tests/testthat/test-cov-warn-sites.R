# The covariance-warning sites (functions that call .covWarn() or .covClip()),
# driven with every label a fit really carries: a label is never what gets
# clipped, so each warning is one whole line of at most .covWarnWidth columns.
# Clipping is for input no label produces (a long condition message), which the
# sites with such input are driven with too.  Sites that need a fit are driven in
# the file of their method; test-cov-warning-gate.R fails for a site no test
# drives.

test_that("the labels come from the code that names them", {
  .fd <- .covFdLabels()
  expect_setequal(
    .fd,
    c(
      c("r", "r+", "|r|", "s", "s+", "|s|"),
      as.vector(outer(c("r", "r+", "|r|"), c("s", "s+", "|s|"), paste, sep = ","))
    )
  )
  expect_identical(
    unname(vapply(.fd, .covFdType, "")),
    ifelse(grepl(",", .fd, fixed = TRUE), "r,s", ifelse(grepl("r", .fd, fixed = TRUE), "r", "s"))
  )
  .labels <- .covRealLabels()
  expect_true(all(
    c("r,s", "r", "s", "analytic", "linFim", "Ha", "|r|,s+", "r+,|s|", "|r|,|s|", "|Ha|", "r,s (full)") %in% .labels
  ))
  expect_false(anyDuplicated(.labels) > 0L)
  expect_false(any(c("", NA_character_) %in% .labels))
  expect_true(all(.covRealLabels(decorated = FALSE, full = FALSE) %in% .labels))
  expect_false(any(vapply(.covRealLabels(full = FALSE), .covIsFull, logical(1))))
  expect_false(any(grepl("|", .covRealLabels(decorated = FALSE), fixed = TRUE)))
  .rep <- .covRealRepaired()
  expect_identical(.covUnrepaired(.rep$label), .rep$request)
  expect_identical(
    .covRealAnalyticKept()[1:3],
    c("r (analytic)", "r+ (analytic)", "|r| (analytic)")
  )
})

test_that(".covInstall() says what it installed instead of the request, whole, for every label", {
  # the grid of label pairs is large
  skip_on_cran()
  .labels <- .covRealLabels()
  .got <- character()
  .want <- character()
  .covSite(".covInstall", {
    for (.have in .labels) {
      for (.asked in .labels) {
        if (.covSameName(.covBaseName(.asked), .covBaseName(.have))) {
          next
        }
        .got <- c(.got, .covSiteWarning(.covInstall(.covSiteEnv(), .covSitePd(), .have, what = .asked)))
        .want <- c(.want, sprintf("\"%s\" installed instead of the requested \"%s\"", .have, .asked))
      }
    }
  })
  expect_gt(length(.want), 1000L)
  expect_identical(.got, .want)
  expect_true(.covSiteWhole(.got))
})

test_that(".covRejectWarn() is one whole line for every label, reason and kept covariance", {
  # the grid of label pairs is large
  skip_on_cran()
  .labels <- .covRealLabels()
  .fd <- .covFdLabels()
  .fdAll <- c(.fd, paste0(.fd, .covFullSuffix))
  .guard <- c("could not be computed", "is not finite", "is not symmetric", "is not positive definite")
  # the guard of .covInstall(): any label rejected, any one kept; the full
  # finite-difference step; the SAEM "Ha" block
  .grid <- rbind(
    expand.grid(
      what = .labels,
      why = .guard,
      kept = c(NA_character_, .labels, .covRealAnalyticKept()),
      stringsAsFactors = FALSE
    ),
    expand.grid(
      what = paste0(c("r", "s", "r,s"), .covFullSuffix),
      why = c("needs a positive-definite R", "needs a positive-definite S"),
      kept = .fdAll,
      stringsAsFactors = FALSE
    ),
    expand.grid(
      what = .saemHaThetaName,
      why = c("has information rows not matching the thetas", "has no estimated mu-referenced theta", .guard),
      kept = c(NA_character_, "linFim", "|linFim|", "linFim (full)"),
      stringsAsFactors = FALSE
    )
  )
  .reject <- function(what, why, kept) {
    .e <- new.env(parent = emptyenv())
    if (!is.na(kept)) {
      .e$cov <- .covSitePd()
      .e$covMethod <- kept
    }
    .covSiteWarning(.covRejectWarn(.e, what, why))
  }
  .got <- .covSite(".covRejectWarn", mapply(.reject, .grid$what, .grid$why, .grid$kept, USE.NAMES = FALSE))
  expect_gt(length(.got), 10000L)
  expect_identical(
    .got,
    sprintf(
      "\"%s\" %s; %s",
      .grid$what,
      .grid$why,
      ifelse(is.na(.grid$kept), "none installed", sprintf("kept \"%s\"", .grid$kept))
    )
  )
  expect_true(.covSiteWhole(.got))
})

test_that(".covRepairWarn() is one whole line for every repaired covariance", {
  .rep <- .covRealRepaired()
  .got <- .covSite(
    ".covRepairWarn",
    mapply(
      function(request, label) .covSiteWarning(.covRepairWarn(request, label)),
      .rep$request,
      .rep$label,
      USE.NAMES = FALSE
    )
  )
  expect_identical(
    .got,
    sprintf("\"%s\" is not positive definite; installed as \"%s\" (sqrtm)", .rep$request, .rep$label)
  )
  expect_true(.covSiteWhole(.got))
})

test_that(".covNotUsableMsg() is one whole line for every walk, request, reason and installed label", {
  .fd <- .covFdLabels()
  .fdAll <- c(.fd, paste0(.fd, .covFullSuffix))
  .tails <- function(labels) c("no covariance", sprintf("installed \"%s\"", labels))
  # FOCEi's finite-difference walk, the after-fit fallback and SAEM's chain, one
  # reason at a time: several reasons joined by "; " can pass .covWarnWidth and are
  # clipped by design (covTried has them in full), as is the analytic reason,
  # which test-covSelect.R pins
  .walks <- list(
    list(
      request = c("r,s", "r", "s", "r,s (full)", "r (full)", "s (full)"),
      why = c(
        "R not computed",
        "R not PD",
        "R failed",
        "S not computed",
        "S not PD",
        "S failed",
        "sandwich doubtful",
        "all variances < 1e-7"
      ),
      tail = .tails(.fdAll)
    ),
    list(
      request = c("r,s", "r", "s", "analytic"),
      why = "could not be computed",
      tail = .tails(c("sa", "imp", "|imp|", "nlme", "linFim"))
    ),
    list(
      request = c("linFim", "linFim (full)", "sa", .saemHaThetaName),
      why = "not usable",
      tail = .tails(c("linFim", "|linFim|", "linFim (full)", "sa", .saemHaThetaName, .saemHaThetaRepairedName))
    )
  )
  .grid <- do.call(
    rbind,
    lapply(.walks, function(w) expand.grid(request = w$request, why = w$why, tail = w$tail, stringsAsFactors = FALSE))
  )
  .got <- .covSite(".covNotUsableMsg", mapply(.covNotUsableMsg, .grid$request, .grid$why, .grid$tail, USE.NAMES = FALSE))
  expect_identical(.got, sprintf("\"%s\" not usable (%s); %s", .grid$request, .grid$why, .grid$tail))
  expect_true(.covSiteWhole(.got))
  expect_gt(length(.got), 500L)
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
        .got <- c(.got, .covSiteWarning(.foceiAnalyticErrWarn(.site)(simpleError(.m))))
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
  expect_true(all(cli::ansi_nchar(.got, type = "width") <= .covWarnWidth))
})

test_that(".foceiWarnConditionalCov() names every finite-difference label whole", {
  .fd <- .covFdLabels()
  .fdAll <- c(.fd, paste0(.fd, .covFullSuffix))
  .conditional <- function(label) {
    .e <- new.env(parent = emptyenv())
    .e$etaObf <- 1
    .e$cov <- diag(2)
    .e$covMethod <- label
    .covSiteWarning(.foceiWarnConditionalCov(.e, list(maxInnerIterations = 0L)))
  }
  .got <- .covSite(".foceiWarnConditionalCov", vapply(.fdAll, .conditional, "", USE.NAMES = FALSE))
  expect_identical(
    .got,
    sprintf("\"%s\" is conditional on the ETAs; setCov() gives the marginal", .fdAll)
  )
  expect_true(.covSiteWhole(.got))
})

test_that(".setCovFd() names the label it installed and the one asked for, whole", {
  .fd <- .covFdLabels()
  .fdAll <- c(.fd, paste0(.fd, .covFullSuffix))
  .cur <- NULL
  local_mocked_bindings(
    .setCovRefit = function(...) list(covMethod = .cur, cov = .covSitePd(), parFixedDf = NULL, parFixed = NULL),
    .covInstall = function(...) TRUE,
    .covCacheDrop = function(...) NULL
  )
  .got <- character()
  .want <- character()
  .covSite(".setCovFd", {
    for (.asked in .fdAll) {
      for (.have in .fdAll) {
        if (.covSameName(.asked, .have) || !identical(.covIsFull(.asked), .covIsFull(.have))) {
          next
        }
        .cur <- .have
        .ctl <- list(covFallback = stats::setNames(list(.fd), .covBaseName(.asked)))
        .got <- c(.got, .covSiteWarning(.setCovFd(NULL, new.env(), .asked, .covBaseName(.asked), control = .ctl)))
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
