# A declined covMethod="analytic" says why.  Each step of the live hook
# (.foceiCalcRanalytic) that gives up records its reason in full (with what the
# failing solve said) in the message and at the head of fit$env$covTried, and the
# step alone, on one line, in the line the C++ hook prints and in $runInfo;
# setCov(fit, "analytic") stops with the full reason.  A path no model reaches on
# its own is forced by replacing the helper behind it.  Every decline path, forced
# one at a time, is enumerated in test-cov-analytic-declines.R.

nmTest({
  .arMod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  # The covariance at the initial estimates is enough to reach every step.
  .arFit <- function(model = .arMod, data = nlmixr2data::theo_sd, maxOuterIterations = 0L, est = "focei", ...) {
    .acc <- new.env(parent = emptyenv())
    .acc$msg <- character(0)
    .out <- utils::capture.output(
      .fit <- withCallingHandlers(
        suppressWarnings(nlmixr2(
          model,
          data,
          est,
          foceiControl(covMethod = "analytic", print = 0L, maxOuterIterations = maxOuterIterations, ...),
          table = tableControl(cwres = FALSE, npde = FALSE)
        )),
        message = function(m) {
          .m <- conditionMessage(m)
          .acc$msg <- c(.acc$msg, sub("\n$", "", .m))
          # nmSuppressMsg() silences the C++ lines when its probe message(" ") is
          # muffled, so let the probe through to keep the printed covariance line
          if (!identical(.m, " \n")) {
            invokeRestart("muffleMessage")
          }
        }
      )
    )
    list(fit = .fit, msg = .acc$msg, out = .out)
  }

  # `reason` is the whole recorded reason; `msg` the message that gave it.
  # the analytic warnings in $runInfo, each of which must stay under 75 characters
  .arWarnings <- function(r) {
    .w <- grep("analytic", r$fit$runInfo, value = TRUE, fixed = TRUE)
    expect_true(all(nchar(.w) < 75L), label = paste(.w, collapse = " | "))
    .w
  }

  # `reason` is the full recorded reason, `short` the step alone, `msg` the message
  .arExpectReason <- function(r, reason, msg, short = reason) {
    expect_false(identical(.covBaseName(r$fit$covMethod), "analytic"))
    expect_identical(r$fit$env$covTried$method[1], "analytic")
    expect_identical(r$fit$env$covTried$outcome[1], reason)
    expect_identical(sum(grepl("covType=\"analytic\": ", r$msg, fixed = TRUE)), 1L)
    expect_true(msg %in% r$msg)
    expect_true(any(grepl(
      paste0(
        "covType=\"analytic\" not used (",
        .foceiAnalyticClip(short),
        "); using the finite-difference \"r,s\" covariance."
      ),
      r$out,
      fixed = TRUE
    )))
    # the warning names the analytic step and then any finite-difference matrix that
    # was not usable either, as covTried does after its first row
    .tried <- r$fit$env$covTried$outcome[-1L]
    .why <- paste(unique(c(short, .tried[.tried != "used"])), collapse = "; ")
    .w <- grep("^\"analytic\" not usable \\(", .arWarnings(r), value = TRUE)
    expect_length(.w, 1L)
    expect_identical(.w, .covNotUsableMsg("analytic", .why, sub("^.*\\); ", "", .w)))
    expect_false(any(grepl("not available for this model", c(r$out, r$fit$runInfo), fixed = TRUE)))
  }

  .arScope <- function(reason) {
    paste0(
      "covType=\"analytic\": ",
      reason,
      " is out of analytic-covariance scope; using the finite-difference covariance instead"
    )
  }

  .arFailed <- function(reason) {
    paste0("covType=\"analytic\": ", reason, "; using the finite-difference covariance instead")
  }

  test_that("an error model the analytic covariance cannot read says so", {
    skip_on_cran()
    local_mocked_bindings(.foceiAnalyticErrFull = function(ui) NULL)
    .r <- .arFit()
    .why <- "a residual error model the analytic covariance cannot read"
    .arExpectReason(.r, .why, .arScope(.why))
  })

  test_that("the reason the error model gave is kept, and said once", {
    skip_on_cran()
    local_mocked_bindings(.foceiAnalyticErrFull = function(ui) .foceiAnalyticFallback("a test error model"))
    .r <- .arFit()
    .arExpectReason(.r, "a test error model", .arScope("a test error model"))
  })

  test_that("a mu-referenced occasion eta says so", {
    skip_on_cran()
    .testSeed(42)
    .sim <- rxode2::rxode2({
      ka <- exp(0.5 + eta.ka)
      cl <- exp(1.0 + eta.cl + iov.cl)
      v <- exp(3.4 + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
    })
    .ev <- rxode2::et(seq(0.25, 20, by = 2))
    .ev <- rxode2::et(.ev, seq(48.25, 68, by = 2))
    .ev <- rxode2::et(.ev, amt = 320, time = 0, cmt = "depot")
    .ev <- as.data.frame(rxode2::et(.ev, amt = 320, time = 48, cmt = "depot"))
    .d <- do.call(
      rbind,
      lapply(seq_len(12), function(id) {
        .dd <- .ev
        .dd$iov.cl <- ifelse(.dd$time < 40, rnorm(1, 0, sqrt(0.08)), rnorm(1, 0, sqrt(0.08)))
        .s <- rxode2::rxSolve(
          .sim,
          c(eta.ka = rnorm(1, 0, sqrt(0.3)), eta.cl = rnorm(1, 0, sqrt(0.1)), eta.v = rnorm(1, 0, sqrt(0.1))),
          .dd,
          returnType = "data.frame"
        )
        .obs <- .s[!is.na(.s$cp) & .s$cp > 0, ]
        data.frame(
          ID = id,
          TIME = .obs$time,
          DV = .obs$cp * exp(rnorm(nrow(.obs), 0, 0.1)),
          AMT = 0,
          EVID = 0,
          occ = ifelse(.obs$time < 40, 1L, 2L)
        )
      })
    )
    .dose <- data.frame(
      ID = rep(1:12, each = 2),
      TIME = rep(c(0, 48), 12),
      DV = 0,
      AMT = 320,
      EVID = 1,
      occ = rep(1:2, 12)
    )
    .d <- rbind(.dose, .d)
    .d <- .d[order(.d$ID, .d$TIME, -.d$EVID), ]
    .iov <- function() {
      ini({
        tka <- 0.5
        tcl <- 1.0
        tv <- 3.4
        eta.ka ~ 0.3
        eta.cl ~ 0.1
        eta.v ~ 0.1
        iov.cl ~ 0.08 | occ
        add.sd <- 0.3
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + iov.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .map <- .foceiEtaThetaMap
    local_mocked_bindings(.foceiEtaThetaMap = function(ui) {
      .m <- .map(ui)
      .m$thetaForEta[grepl("^rx\\.iov\\.cl\\.", .m$etaNames)] <- "tcl"
      .m
    })
    .r <- .arFit(.iov, .d, iovXform = "sd")
    .arExpectReason(.r, "a mu-referenced occasion (IOV) eta", .arScope("a mu-referenced occasion (IOV) eta"))
  })

  test_that("covariance parameters missing from the assembly are named", {
    skip_on_cran()
    .dirs <- .foceiAnalyticDirections
    local_mocked_bindings(.foceiAnalyticDirections = function(...) {
      .d <- .dirs(...)
      .d$sgName <- character(0)
      .d
    })
    .r <- .arFit()
    .why <- "covariance parameters the analytic assembly does not carry (add.sd)"
    .arExpectReason(.r, .why, .arScope(.why))
  })

  test_that("an augmented model that will not build says so", {
    skip_on_cran()
    local_mocked_bindings(.foceiAnalyticAugModelDirs = function(ui, dirs, order = 2L) NULL)
    .r <- .arFit()
    .why <- "an augmented model that would not build"
    .arExpectReason(.r, .why, .arFailed(.why))
  })

  test_that("a failed subject solve names the subject and what the solve said", {
    skip_on_cran()
    # the batched solve fails first; its failure is recovered, so only the
    # per-subject one is reported
    local_mocked_bindings(
      .foceiAnalyticSolveAllFD3 = function(...) .foceiAnalyticNote("a batched failure"),
      .foceiAnalyticSolveSubjectFD3 = function(...) .foceiAnalyticNote("a test solve failure")
    )
    .r <- .arFit()
    .why <- "a failed sensitivity solve for subject 1 (a test solve failure)"
    .arExpectReason(.r, .why, .arFailed(.why), "a failed sensitivity solve for subject 1")
  })

  test_that("residual-error parameters the assembly did not order say so", {
    skip_on_cran()
    local_mocked_bindings(.foceiAnalyticSgNames = function(dir, ef, general) NULL)
    .r <- .arFit()
    .why <- "residual-error parameters the assembly did not order"
    .arExpectReason(.r, .why, .arFailed(.why))
  })

  test_that("an observed information that will not invert gives what solve() said", {
    skip_on_cran()
    local_mocked_bindings(.foceiAnalyticSolve = function(m) .foceiAnalyticNote("a test singular matrix"))
    .r <- .arFit()
    .why <- "an observed information that would not invert (a test singular matrix)"
    .arExpectReason(.r, .why, .arFailed(.why), "an observed information that would not invert")
  })

  test_that("a theta block that will not invert says so", {
    skip_on_cran()
    .calls <- new.env(parent = emptyenv())
    .calls$n <- 0L
    local_mocked_bindings(.foceiAnalyticSolve = function(m) {
      .calls$n <- .calls$n + 1L
      if (.calls$n == 1L) solve(m) else .foceiAnalyticNote("a test singular block")
    })
    .r <- .arFit()
    expect_identical(.calls$n, 2L)
    .why <- "a theta block of the covariance that would not invert (a test singular block)"
    .arExpectReason(.r, .why, .arFailed(.why), "a theta block of the covariance that would not invert")
  })

  test_that("an error in the assembly is the reason, and still a warning", {
    skip_on_cran()
    local_mocked_bindings(.foceiAnalyticDirections = function(...) stop("a test error"))
    .r <- .arFit()
    expect_identical(.r$fit$env$covTried$outcome[1], "an error (a test error)")
    expect_true("analytic err1: a test error" %in% .arWarnings(.r))
    expect_true(any(grepl("covType=\"analytic\" not used (an error)", .r$out, fixed = TRUE)))
    expect_true(any(startsWith(.arWarnings(.r), "\"analytic\" not usable (an error); installed \"")))
    # a long, multi-line error is one line, and clipped in the warning
    local_mocked_bindings(.foceiAnalyticDirections = function(...) {
      stop(paste(rep("a long test error", 6), collapse = "\n"))
    })
    .r <- .arFit()
    expect_identical(
      .r$fit$env$covTried$outcome[1],
      paste0("an error (", paste(rep("a long test error", 6), collapse = " "), ")")
    )
    expect_true(
      paste0("analytic err1: ", .foceiAnalyticClip(paste(rep("a long test error", 6), collapse = " "))) %in%
        .arWarnings(.r)
    )
  })

  test_that("with no fallback listed, the failed covariance still says why", {
    skip_on_cran()
    local_mocked_bindings(.foceiAnalyticSgNames = function(dir, ef, general) NULL)
    .r <- .arFit(covFallback = list())
    .why <- "residual-error parameters the assembly did not order"
    expect_identical(.r$fit$covMethod, "failed")
    expect_identical(.r$fit$env$covTried, data.frame(method = "analytic", outcome = .why))
    expect_true(any(grepl(
      paste0("covType=\"analytic\" not used (", .why, "), and covFallback lists no fallback for it."),
      .r$out,
      fixed = TRUE
    )))
    expect_true(.covNotUsableMsg("analytic", .why, "covFallback lists no fallback") %in% .arWarnings(.r))
  })

  test_that("an analytic covariance that is used records no decline", {
    skip_on_cran()
    # converged: at the initial estimates the full matrix is not positive definite
    .r <- .arFit(maxOuterIterations = 5000L)
    expect_identical(.covBaseName(.r$fit$covMethod), "analytic")
    expect_identical(.r$fit$env$covTried, data.frame(method = "r", outcome = "used"))
    expect_false(exists(".analyticReason", envir = .r$fit$env, inherits = FALSE))
  })

  # add + prop error: the add/prop assembly, whose FOCE freezes R at an eta = 0 solve
  .arModAddProp <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.5
      prop.sd <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) + prop(prop.sd)
    })
  }

  # combined1: the general (f,R) assembly
  .arModC1 <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.5
      prop.sd <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) + prop(prop.sd) + combined1()
    })
  }

  # IDs 101..112, so a reason naming a subject shows the data's ID, not the internal code
  .arData101 <- nlmixr2data::theo_sd
  .arData101$ID <- .arData101$ID + 100L

  test_that("a failed eta = 0 population solve says so, for both FOCE assemblies", {
    skip_on_cran()
    local_mocked_bindings(.foceiAnalyticSolveAll = function(...) NULL)
    .why <- "a failed eta=0 population solve"
    for (.m in list(.arModAddProp, .arModC1)) {
      .arExpectReason(.arFit(.m, est = "foce"), .why, .arFailed(.why))
    }
  })

  test_that("the population solve names the subject and what went wrong", {
    skip_on_cran()
    .rx <- .foceiAnalyticRxSolve
    # every row of the first subject gone from the population solve, or its first row
    local_mocked_bindings(.foceiAnalyticRxSolve = function(...) {
      .d <- .rx(...)
      if (is.null(.d) || is.null(.d$id)) .d else .d[as.character(.d$id) != "1", , drop = FALSE]
    })
    .why <- "a failed eta=0 population solve (no solved rows for subject 101)"
    .arExpectReason(
      .arFit(.arModAddProp, .arData101, est = "foce"),
      .why,
      .arFailed(.why),
      "a failed eta=0 population solve"
    )
    local_mocked_bindings(.foceiAnalyticRxSolve = function(...) {
      .d <- .rx(...)
      if (is.null(.d)) .d else .d[-1L, , drop = FALSE]
    })
    .why <- "a failed eta=0 population solve (subject 101: 10 observation rows solved for 11 observations)"
    .arExpectReason(
      .arFit(.arModAddProp, .arData101, est = "foce"),
      .why,
      .arFailed(.why),
      "a failed eta=0 population solve"
    )
  })

  test_that("a subject solve whose rows do not line up says how, for each way", {
    skip_on_cran()
    .rx <- .foceiAnalyticRxSolve
    .alter <- list(
      "the augmented solve returned no observation rows" = function(d) d[0L, , drop = FALSE],
      "10 observation rows solved for 11 observations" = function(d) d[-1L, , drop = FALSE],
      "solved observation times differ from the data's" = function(d) {
        d$time[1L] <- d$time[1L] + 1e-3
        d
      }
    )
    for (.n in names(.alter)) {
      local_mocked_bindings(.foceiAnalyticRxSolve = function(...) {
        .d <- .rx(...)
        if (is.null(.d)) .d else .alter[[.n]](.d)
      })
      .short <- "a failed sensitivity solve for subject 101"
      .why <- paste0(.short, " (", .n, ")")
      .arExpectReason(.arFit(data = .arData101), .why, .arFailed(.why), .short)
    }
  })

  test_that("a failed or non-finite observed-information kernel says so", {
    skip_on_cran()
    for (.k in list(list(est = "focei", fn = "foceiRAllFR_"), list(est = "foce", fn = "foceiRAllFoceFR_"))) {
      local({
        .m <- list(function(...) stop("a test kernel failure"))
        names(.m) <- .k$fn
        do.call(local_mocked_bindings, .m)
        .short <- "a failed observed-information kernel"
        .why <- paste0(.short, " (a test kernel failure)")
        .arExpectReason(.arFit(.arModC1, est = .k$est), .why, .arFailed(.why), .short)
      })
      local({
        .m <- list(function(...) matrix(NaN, 2L, 2L))
        names(.m) <- .k$fn
        do.call(local_mocked_bindings, .m)
        .why <- "a non-finite observed information"
        .arExpectReason(.arFit(.arModC1, est = .k$est), .why, .arFailed(.why))
      })
    }
  })

  test_that("a subject's observed information that fails or is not finite names the subject", {
    skip_on_cran()
    local({
      local_mocked_bindings(.foceiAnalyticSubjectR = function(...) NULL)
      .why <- "an observed information that would not assemble for subject 101"
      .arExpectReason(.arFit(data = .arData101), .why, .arFailed(.why))
    })
    local({
      local_mocked_bindings(.foceiAnalyticSubjectR = function(...) matrix(NaN, 1L, 1L))
      .why <- "a non-finite observed information for subject 101"
      .arExpectReason(.arFit(data = .arData101), .why, .arFailed(.why))
    })
  })

  test_that("a failed quadrature-node solve names the subject and what the solve said", {
    skip_on_cran()
    # the node solves are the only subject solves left once the batched solve works
    local_mocked_bindings(.foceiAnalyticSolveFA = function(...) .foceiAnalyticNote("a test node failure"))
    .short <- "an observed information that would not assemble for subject 101"
    .why <- paste0(.short, " (a failed quadrature-node solve: a test node failure)")
    .arExpectReason(.arFit(data = .arData101, nAGQ = 3L), .why, .arFailed(.why), .short)
  })

  test_that("setCov(fit, \"analytic\") stops with the reason", {
    skip_on_cran()
    .fit <- suppressMessages(suppressWarnings(nlmixr2(
      .arMod,
      nlmixr2data::theo_sd,
      "focei",
      foceiControl(covMethod = "r,s", print = 0L),
      table = tableControl(cwres = FALSE, npde = FALSE)
    )))
    .before <- .fit$covMethod
    .fail <- function(why) {
      paste0(
        "covMethod=\"analytic\" could not be computed for this fit (",
        why,
        "); the covariance is left unchanged"
      )
    }
    local({
      local_mocked_bindings(.foceiAnalyticSolve = function(m) .foceiAnalyticNote("a test singular matrix"))
      expect_error(
        suppressMessages(setCov(.fit, "analytic")),
        .fail("an observed information that would not invert (a test singular matrix)"),
        fixed = TRUE
      )
    })
    local({
      local_mocked_bindings(.foceiAnalyticErrFull = function(ui) NULL)
      expect_error(
        suppressMessages(setCov(.fit, "analytic")),
        .fail("a residual error model the analytic covariance cannot read"),
        fixed = TRUE
      )
    })
    local({
      local_mocked_bindings(.foceiAnalyticDirections = function(...) stop("a test error"))
      expect_error(
        suppressMessages(setCov(.fit, "analytic")),
        .fail("an error (a test error)"),
        fixed = TRUE
      )
    })
    expect_identical(.fit$covMethod, .before)
  })

  test_that("an augmented model that will not build gives the build's error", {
    skip_on_cran()
    # the build is cached per model, so empty the cache to make this one build again
    rm(list = ls(.foceiAnalyticAugCache, all.names = TRUE), envir = .foceiAnalyticAugCache)
    local_mocked_bindings(.foceiAnalyticCols = function(...) stop("a test build error"))
    .short <- "an augmented model that would not build"
    .why <- paste0(.short, " (a test build error)")
    .arExpectReason(.arFit(), .why, .arFailed(.why), .short)
  })
})
