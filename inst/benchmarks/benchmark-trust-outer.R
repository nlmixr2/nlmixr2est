#!/usr/bin/env Rscript
# Benchmark est="trust" (the C++-resident RcppTrust outer optimizer for the
# population theta vector, R/trust.R + src/nlm.cpp's nlmTrustFit()) against
# est="bobyqa" (the nlm-family default derivative-free optimizer) across every
# population-only ("nlm style") model fit in this package's OWN test suite --
# not the nlmixr2 vignettes (those are FOCEi/mixed-effects models; see
# benchmark-trust-inner.R for that corpus). This is a DIFFERENT axis from
# benchmark-trust-inner.R: that script compares foceiControl(innerOpt=) (the
# per-subject ETA optimizer inside FOCEi); this one compares which OUTER
# optimizer drives the population THETA vector for nlm-family fits.
#
# Corpus: every tests/testthat/*.R file found (as of this writing) to use one
# of the nlm-family est= methods (nlm/nlminb/bobyqa/newuoa/uobyqa/n1qn1/
# lbfgsb3c/optim) --
#   test-cens-dist-warn.R, test-lincmt-ode-fit.R, test-matexp.R,
#   test-nlm-cens-t.R, test-nlm-lik-contrib.R, test-nlm.R,
#   test-nlmsetup-fresh-rx.R, test-optexpression-ebe.R, test-optim.R,
#   test-splitbolus-interp.R
# (test-trust.R itself is excluded -- it is this feature's own unit-test
# file, with a single small smoke-test model already covered incidentally by
# every other model in this corpus.)
#
# Like benchmark-trust-inner.R's corpus (b)/(c), the real test-file source is
# executed with nlmixr2()/nlmixr() intercepted: whenever a call resolves to a
# nlm-family fit (any est other than "trust" itself -- none of the source
# files predate/reference "trust", but the check is symmetric on principle),
# it is re-run once with est="bobyqa" and once with est="trust" (fresh
# default bobyqaControl()/trustControl(), NOT the original call's own control
# object -- nlm-family control objects are per-method R functions with
# incompatible field sets, so there is no generic way to "carry over" one
# method's control fields into another's the way foceiControl(innerOpt=)'s
# single shared control object lets benchmark-trust-inner.R do; a clean
# default-vs-default comparison is the fair, well-defined thing to measure
# instead). This also naturally discovers whether every claimed
# "population-only" model in the corpus really is one -- if a model turns out
# to have an eta the original test's est= happened to silently ignore, the
# re-fit either succeeds (population-only in practice) or errors (logged, not
# silently dropped, per this repo's own benchmark-harness convention).
#
# Usage:
#   Rscript inst/benchmarks/benchmark-trust-outer.R [outDir] [testDir]
# outDir defaults to inst/benchmarks/results; testDir defaults to
# tests/testthat.

suppressPackageStartupMessages({
  library(nlmixr2est)
})

.bargs <- commandArgs(trailingOnly = TRUE)
.outDir <- if (length(.bargs) >= 1) .bargs[1] else "inst/benchmarks/results"
.testDir <- if (length(.bargs) >= 2) .bargs[2] else "tests/testthat"
dir.create(.outDir, showWarnings = FALSE, recursive = TRUE)
.csvPath <- file.path(.outDir, "trust-outer-benchmark.csv")
.logPath <- file.path(.outDir, "trust-outer-benchmark.log")

.logCon <- file(.logPath, open = "a")
.blog <- function(fmt, ...) {
  msg <- sprintf(fmt, ...)
  cat(msg, "\n")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg, "\n", file = .logCon)
  flush(.logCon)
}
.blog("=== benchmark-trust-outer.R starting ===")

# ---------------------------------------------------------------------------
# Result accumulation (incremental CSV write -- a crash mid-sweep still
# leaves every completed row on disk).
# ---------------------------------------------------------------------------
.hessianMethods <- c("fd", "bfgs", "sr1", "bofill")

.benchResults <- list()
.benchRecord <- function(source, model, bq, trList) {
  row <- data.frame(
    source = source, model = model,
    bobyqa_ok = isTRUE(bq$ok), bobyqa_time = as.numeric(bq$time), bobyqa_objf = as.numeric(bq$objf),
    bobyqa_error = ifelse(is.null(bq$error), NA_character_, bq$error),
    stringsAsFactors = FALSE
  )
  logParts <- sprintf("bobyqa ok=%s time=%.2fs objf=%s", bq$ok, bq$time, format(bq$objf))
  for (hm in .hessianMethods) {
    tr <- trList[[hm]]
    maxParDiff <- NA_real_
    if (isTRUE(bq$ok) && isTRUE(tr$ok) &&
          !is.null(bq$theta) && !is.null(tr$theta)) {
      common <- intersect(names(bq$theta), names(tr$theta))
      if (length(common)) maxParDiff <- max(abs(bq$theta[common] - tr$theta[common]))
    }
    row[[paste0("trust_", hm, "_ok")]] <- isTRUE(tr$ok)
    row[[paste0("trust_", hm, "_time")]] <- as.numeric(tr$time)
    row[[paste0("trust_", hm, "_objf")]] <- as.numeric(tr$objf)
    row[[paste0("trust_", hm, "_error")]] <- ifelse(is.null(tr$error), NA_character_, tr$error)
    row[[paste0("trust_", hm, "_dObjf")]] <- as.numeric(tr$objf) - as.numeric(bq$objf)
    row[[paste0("trust_", hm, "_speedup")]] <- as.numeric(bq$time) / as.numeric(tr$time)
    row[[paste0("trust_", hm, "_maxParDiff")]] <- maxParDiff
    logParts <- c(logParts, sprintf("trust[%s] ok=%s time=%.2fs objf=%s dObjf=%s maxParDiff=%s",
                                    hm, tr$ok, tr$time, format(tr$objf),
                                    format(row[[paste0("trust_", hm, "_dObjf")]]), format(maxParDiff)))
  }
  .benchResults[[length(.benchResults) + 1]] <<- row
  write.csv(do.call(rbind, .benchResults), .csvPath, row.names = FALSE)
  .blog("[%s] %s: %s", source, model, paste(logParts, collapse = " | "))
}

.benchLogSkip <- function(source, model, reason) {
  .blog("[%s] %s: SKIPPED -- %s", source, model, reason)
  row <- data.frame(
    source = source, model = model,
    bobyqa_ok = NA, bobyqa_time = NA_real_, bobyqa_objf = NA_real_, bobyqa_error = reason,
    stringsAsFactors = FALSE
  )
  for (hm in .hessianMethods) {
    row[[paste0("trust_", hm, "_ok")]] <- NA
    row[[paste0("trust_", hm, "_time")]] <- NA_real_
    row[[paste0("trust_", hm, "_objf")]] <- NA_real_
    row[[paste0("trust_", hm, "_error")]] <- reason
    row[[paste0("trust_", hm, "_dObjf")]] <- NA_real_
    row[[paste0("trust_", hm, "_speedup")]] <- NA_real_
    row[[paste0("trust_", hm, "_maxParDiff")]] <- NA_real_
  }
  .benchResults[[length(.benchResults) + 1]] <<- row
  write.csv(do.call(rbind, .benchResults), .csvPath, row.names = FALSE)
}

# ---------------------------------------------------------------------------
# Core timed-fit helper: refits (object, data) fresh under the given est=,
# with a plain default control (print=0, calcTables=FALSE) for that method.
# ---------------------------------------------------------------------------
.benchFitOnce <- function(object, data, est, table, dots, save, envir, hessianMethod = "fd") {
  ctl <- switch(est,
    bobyqa = nlmixr2est::bobyqaControl(print = 0L, calcTables = FALSE),
    trust  = nlmixr2est::trustControl(print = 0L, calcTables = FALSE, hessianMethod = hessianMethod),
    stop("unsupported est for .benchFitOnce: ", est)
  )
  t0 <- proc.time()["elapsed"]
  fit <- tryCatch(
    suppressWarnings(suppressMessages(
      do.call(nlmixr2est::nlmixr2,
              c(list(object = object, data = data, est = est, control = ctl,
                     table = table, save = save, envir = envir), dots))
    )),
    error = function(e) e
  )
  t1 <- proc.time()["elapsed"]
  if (inherits(fit, "error")) {
    return(list(ok = FALSE, time = unname(t1 - t0), objf = NA_real_, theta = NULL,
                error = conditionMessage(fit)))
  }
  list(ok = TRUE, time = unname(t1 - t0),
       objf = tryCatch(as.numeric(fit$objective), error = function(e) NA_real_),
       theta = tryCatch(setNames(as.numeric(fit$theta), names(fit$theta)), error = function(e) NULL),
       error = NA_character_)
}

# ---------------------------------------------------------------------------
# Interception wrapper installed as nlmixr2()/nlmixr() while sourcing corpus
# files. Runs the ORIGINAL call once (so downstream code in the sourced file
# still sees a normal result and its own assertions still run); if the
# resolved fit is a nlm-family method OTHER than "trust" itself, additionally
# times a fresh bobyqa run (skipped when the original already used bobyqa)
# and a fresh trust run, and records both.
# ---------------------------------------------------------------------------
.nlmFamilyEst <- c("nlm", "nlminb", "bobyqa", "newuoa", "uobyqa", "n1qn1", "lbfgsb3c", "optim")

.benchSourceTag <- "unknown"
.benchModelCounter <- 0L
.nextBenchModelLabel <- function() {
  .benchModelCounter <<- .benchModelCounter + 1L
  sprintf("%s#%d", .benchSourceTag, .benchModelCounter)
}

.origNlmixr2 <- nlmixr2est::nlmixr2

.benchNlmixr2Wrap <- function(object, data, est = NULL, control = list(),
                               table = nlmixr2est::tableControl(), ..., save = NULL,
                               envir = parent.frame()) {
  dots <- list(...)
  t0 <- proc.time()["elapsed"]
  fit <- tryCatch(
    do.call(.origNlmixr2,
            c(list(object = object, data = data, est = est, control = control,
                   table = table, save = save, envir = envir), dots)),
    error = function(e) e
  )
  t1 <- proc.time()["elapsed"]
  if (inherits(fit, "error")) return(fit) # let the sourced file's own error handling see it
  resolvedEst <- tryCatch(fit$est, error = function(e) NA_character_)
  isNlmFamily <- isTRUE(resolvedEst %in% .nlmFamilyEst)
  if (isNlmFamily) {
    label <- .nextBenchModelLabel()
    if (identical(resolvedEst, "bobyqa")) {
      bq <- list(ok = TRUE, time = unname(t1 - t0),
                 objf = tryCatch(as.numeric(fit$objective), error = function(e) NA_real_),
                 theta = tryCatch(setNames(as.numeric(fit$theta), names(fit$theta)), error = function(e) NULL),
                 error = NA_character_)
    } else {
      bq <- .benchFitOnce(object, data, "bobyqa", table, dots, save, envir)
    }
    trList <- stats::setNames(
      lapply(.hessianMethods, function(hm) {
        .benchFitOnce(object, data, "trust", table, dots, save, envir, hessianMethod = hm)
      }),
      .hessianMethods
    )
    .benchRecord(.benchSourceTag, label, bq, trList)
  }
  fit
}

.withBenchOverrides <- function(envir) {
  assign("nlmixr2", .benchNlmixr2Wrap, envir = envir)
  assign("nlmixr", .benchNlmixr2Wrap, envir = envir)
  # Neutralize testthat's skip/report machinery -- run the block, don't assert.
  # A genuine error inside the block (e.g. an undefined helper) is LOGGED, not
  # silently discarded -- a swallowed error here means zero models were
  # intercepted, indistinguishable in the CSV from "this file legitimately has
  # no nlm-family fit" unless it is surfaced.
  assign("test_that", function(desc, code) {
    tryCatch(force(code), error = function(e) {
      .blog("[%s] test_that(%s) errored: %s", .benchSourceTag, desc, conditionMessage(e))
    })
  }, envir = envir)
  assign("nmTest", function(code) {
    tryCatch(force(code), error = function(e) {
      .blog("[%s] nmTest block errored: %s", .benchSourceTag, conditionMessage(e))
    })
  }, envir = envir)
  assign("skip_on_cran", function(...) invisible(NULL), envir = envir)
  assign("skip_if_not_installed", function(...) invisible(NULL), envir = envir)
  assign("skip_on_ci", function(...) invisible(NULL), envir = envir)
  assign("skip", function(...) invisible(NULL), envir = envir)
}

# Helper files the test suite's own test_that blocks rely on (e.g. .nlmixr(),
# the suppressWarnings/suppressMessages wrapper every corpus file calls
# instead of nlmixr()/nlmixr2() directly -- .benchNlmixr2Wrap intercepts the
# UNDERLYING call .nlmixr() makes, so sourcing these first is what lets the
# interception reach the corpus files at all). Matches this repo's own
# "source all helpers EXCEPT helper-zzz-fits.R" convention (that file drives
# the (irrelevant, slow) FOCEi fixture cache, not needed here).
.sourceHelpers <- function(envir) {
  helperDir <- "tests/testthat"
  helpers <- sort(list.files(helperDir, pattern = "^helper-.*\\.R$", full.names = TRUE))
  helpers <- helpers[basename(helpers) != "helper-zzz-fits.R"]
  for (h in helpers) {
    tryCatch(source(h, local = envir, echo = FALSE),
             error = function(e) .blog("helper source() failed for %s: %s", h, conditionMessage(e)))
  }
}

.sourceCorpusFile <- function(path, sourceTag) {
  .benchSourceTag <<- sourceTag
  .benchModelCounter <<- 0L
  env <- new.env(parent = globalenv())
  .sourceHelpers(env)
  .withBenchOverrides(env)
  ok <- tryCatch({
    source(path, local = env, echo = FALSE)
    TRUE
  }, error = function(e) {
    .blog("[%s] source() failed: %s", sourceTag, conditionMessage(e))
    FALSE
  })
  if (ok && .benchModelCounter == 0L) {
    .benchLogSkip(sourceTag, "(whole file)", "sourced cleanly but no nlm-family fit was intercepted")
  } else if (!ok) {
    .benchLogSkip(sourceTag, "(whole file)", "source() failed -- see log")
  }
  invisible(ok)
}

# ---------------------------------------------------------------------------
# The nlm-family-using test-suite corpus.
# ---------------------------------------------------------------------------
.runNlmFamilyTests <- function(testDir) {
  if (!dir.exists(testDir)) {
    .blog("test directory not found: %s -- aborting", testDir)
    return(invisible(NULL))
  }
  files <- c("test-cens-dist-warn.R", "test-lincmt-ode-fit.R", "test-matexp.R",
             "test-nlm-cens-t.R", "test-nlm-lik-contrib.R", "test-nlm.R",
             "test-nlmsetup-fresh-rx.R", "test-optexpression-ebe.R",
             "test-optim.R", "test-splitbolus-interp.R")
  files <- file.path(testDir, files)
  missing <- files[!file.exists(files)]
  for (m in missing) .blog("corpus file not found, skipping: %s", m)
  files <- files[file.exists(files)]
  .blog("found %d/%d corpus files", length(files), length(files) + length(missing))
  for (f in files) {
    tag <- paste0("test:", sub("\\.R$", "", basename(f)))
    .sourceCorpusFile(f, tag)
  }
}

# ---------------------------------------------------------------------------
# Run the corpus, then write a markdown summary.
# ---------------------------------------------------------------------------
.blog("--- nlm-family test-suite corpus ---")
.runNlmFamilyTests(.testDir)

if (length(.benchResults)) {
  df <- do.call(rbind, .benchResults)
  write.csv(df, .csvPath, row.names = FALSE)
  mdPath <- file.path(.outDir, "trust-outer-benchmark.md")
  header <- paste0("| source | model | bobyqa ok | bobyqa time | bobyqa objf | ",
                    paste(sapply(.hessianMethods, function(hm)
                      sprintf("trust[%s] ok | trust[%s] time | trust[%s] objf | trust[%s] dObjf | trust[%s] speedup | trust[%s] maxParDiff", hm, hm, hm, hm, hm, hm)),
                      collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", 5 + 6 * length(.hessianMethods)), collapse = "|"), "|")
  lines <- c(
    "# est=\"trust\" hessianMethod= vs est=\"bobyqa\" benchmark",
    "",
    sprintf("Total corpus entries: %d", nrow(df)),
    ""
  )
  for (hm in .hessianMethods) {
    okCol <- paste0("trust_", hm, "_ok")
    ok <- df[!is.na(df$bobyqa_ok) & df$bobyqa_ok & !is.na(df[[okCol]]) & df[[okCol]], , drop = FALSE]
    dObjfCol <- ok[[paste0("trust_", hm, "_dObjf")]]
    speedupCol <- ok[[paste0("trust_", hm, "_speedup")]]
    maxParDiffCol <- ok[[paste0("trust_", hm, "_maxParDiff")]]
    lines <- c(lines,
      sprintf("## hessianMethod=\"%s\" (both converged: %d/%d)", hm, nrow(ok), nrow(df)),
      "")
    if (nrow(ok)) {
      lines <- c(lines,
        sprintf("- Median speedup (bobyqa time / trust time): %.2fx", stats::median(speedupCol, na.rm = TRUE)),
        sprintf("- Median |objf diff|: %.4g", stats::median(abs(dObjfCol), na.rm = TRUE)),
        sprintf("- Median max |param diff|: %.4g", stats::median(maxParDiffCol, na.rm = TRUE)),
        sprintf("- Entries with |objf diff| > 1 (likely a different local optimum, not just numeric noise): %d",
                sum(abs(dObjfCol) > 1, na.rm = TRUE)),
        "")
    }
  }
  lines <- c(lines, header, sep)
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    cells <- c(r$source, r$model, r$bobyqa_ok, format(r$bobyqa_time), format(r$bobyqa_objf))
    for (hm in .hessianMethods) {
      cells <- c(cells,
        r[[paste0("trust_", hm, "_ok")]],
        format(r[[paste0("trust_", hm, "_time")]]),
        format(r[[paste0("trust_", hm, "_objf")]]),
        format(r[[paste0("trust_", hm, "_dObjf")]]),
        format(r[[paste0("trust_", hm, "_speedup")]]),
        format(r[[paste0("trust_", hm, "_maxParDiff")]]))
    }
    lines <- c(lines, paste0("| ", paste(cells, collapse = " | "), " |"))
  }
  writeLines(lines, mdPath)
  .blog("wrote %s and %s", .csvPath, mdPath)
} else {
  .blog("no results recorded at all -- something is wrong with the harness or the corpus")
}
close(.logCon)
