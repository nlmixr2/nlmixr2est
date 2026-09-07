#!/usr/bin/env Rscript
# Benchmark innerOpt="trust" (RcppTrust) against innerOpt="n1qn1" -- the two
# optimizers innerOpt="auto" chooses between -- across three corpora:
#   (a) the 5 outer-optimizer benchmark models from discussion #924
#       (https://github.com/nlmixr2/nlmixr2est/discussions/924)
#   (b) every fitting nlmixr2 package vignette (~/src/nlmixr2/vignettes)
#   (c) every non-wang2007 tests/testthat/test-focei-*.R model
#
# For (b)/(c) the real vignette/test source is executed with nlmixr2()/nlmixr()
# intercepted: whenever a call resolves to a "focei" fit, it is re-run once with
# innerOpt="n1qn1" and once with innerOpt="trust" (same ui/data/other control),
# and both runs are timed and compared. This avoids hand-copying every model --
# the corpus is exercised as authored.
#
# Usage:
#   Rscript inst/benchmarks/benchmark-trust-inner.R [outDir] [vignetteDir]
# outDir defaults to inst/benchmarks/results; vignetteDir defaults to
# ~/src/nlmixr2/vignettes (override for a different checkout location).
# Any corpus entry that errors, fails to fit, or is skipped is logged, not
# silently dropped -- see the .log file next to the CSV.

suppressPackageStartupMessages({
  library(nlmixr2est)
})

.bargs <- commandArgs(trailingOnly = TRUE)
.outDir <- if (length(.bargs) >= 1) .bargs[1] else "inst/benchmarks/results"
.vigDir <- if (length(.bargs) >= 2) .bargs[2] else path.expand("~/src/nlmixr2/vignettes")
dir.create(.outDir, showWarnings = FALSE, recursive = TRUE)
.csvPath <- file.path(.outDir, "trust-inner-benchmark.csv")
.logPath <- file.path(.outDir, "trust-inner-benchmark.log")

.logCon <- file(.logPath, open = "a")
.blog <- function(fmt, ...) {
  msg <- sprintf(fmt, ...)
  cat(msg, "\n")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg, "\n", file = .logCon)
  flush(.logCon)
}
.blog("=== benchmark-trust-inner.R starting ===")

# ---------------------------------------------------------------------------
# Result accumulation (incremental CSV write -- a crash mid-sweep still leaves
# every completed row on disk).
# ---------------------------------------------------------------------------
.benchResults <- list()
.benchRecord <- function(source, model, n1, tr) {
  maxParDiff <- NA_real_
  if (isTRUE(n1$ok) && isTRUE(tr$ok) &&
        !is.null(n1$theta) && !is.null(tr$theta)) {
    common <- intersect(names(n1$theta), names(tr$theta))
    if (length(common)) maxParDiff <- max(abs(n1$theta[common] - tr$theta[common]))
  }
  row <- data.frame(
    source = source, model = model,
    n1qn1_ok = isTRUE(n1$ok), n1qn1_time = as.numeric(n1$time), n1qn1_objf = as.numeric(n1$objf),
    n1qn1_error = ifelse(is.null(n1$error), NA_character_, n1$error),
    trust_ok = isTRUE(tr$ok), trust_time = as.numeric(tr$time), trust_objf = as.numeric(tr$objf),
    trust_error = ifelse(is.null(tr$error), NA_character_, tr$error),
    dObjf = as.numeric(tr$objf) - as.numeric(n1$objf),
    speedup = as.numeric(n1$time) / as.numeric(tr$time),
    maxParDiff = maxParDiff,
    stringsAsFactors = FALSE
  )
  .benchResults[[length(.benchResults) + 1]] <<- row
  write.csv(do.call(rbind, .benchResults), .csvPath, row.names = FALSE)
  .blog("[%s] %s: n1qn1 ok=%s time=%.2fs objf=%s | trust ok=%s time=%.2fs objf=%s | dObjf=%s maxParDiff=%s",
        source, model, n1$ok, n1$time, format(n1$objf), tr$ok, tr$time, format(tr$objf),
        format(row$dObjf), format(maxParDiff))
}

.benchLogSkip <- function(source, model, reason) {
  .blog("[%s] %s: SKIPPED -- %s", source, model, reason)
  .benchResults[[length(.benchResults) + 1]] <<- data.frame(
    source = source, model = model,
    n1qn1_ok = NA, n1qn1_time = NA_real_, n1qn1_objf = NA_real_, n1qn1_error = reason,
    trust_ok = NA, trust_time = NA_real_, trust_objf = NA_real_, trust_error = reason,
    dObjf = NA_real_, speedup = NA_real_, maxParDiff = NA_real_,
    stringsAsFactors = FALSE
  )
  write.csv(do.call(rbind, .benchResults), .csvPath, row.names = FALSE)
}

# ---------------------------------------------------------------------------
# Core timed-fit helper: runs `est="focei"` at a forced innerOpt, on the
# already-resolved ui/data, keeping every other control field the caller used.
# ---------------------------------------------------------------------------
.benchFitOnce <- function(object, data, est, baseControl, table, dots, save, envir, innerOpt) {
  ctlList <- tryCatch(as.list(baseControl), error = function(e) list())
  ctlList$innerOpt <- innerOpt
  ctlList$calcTables <- FALSE
  ctl <- tryCatch(do.call(nlmixr2est::foceiControl, ctlList), error = function(e) NULL)
  if (is.null(ctl)) {
    return(list(ok = FALSE, time = NA_real_, objf = NA_real_, theta = NULL,
                error = "could not rebuild foceiControl for this innerOpt"))
  }
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
       objf = tryCatch(as.numeric(fit$objf), error = function(e) NA_real_),
       theta = tryCatch(setNames(as.numeric(fit$theta), names(fit$theta)), error = function(e) NULL),
       error = NA_character_)
}

# ---------------------------------------------------------------------------
# Interception wrapper installed as nlmixr2()/nlmixr() while sourcing corpus
# (b)/(c) files. Runs the ORIGINAL call once; if it produced a "focei" fit,
# additionally times an explicit n1qn1 baseline (skipped when the original
# already used n1qn1 -- the common case) and a trust run, and records both.
# Always returns the original fit so downstream code in the sourced file sees
# a normal result.
# ---------------------------------------------------------------------------
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
  isFocei <- tryCatch(identical(fit$est, "focei"), error = function(e) FALSE)
  if (isTRUE(isFocei)) {
    label <- .nextBenchModelLabel()
    baseCtl <- tryCatch(fit$control, error = function(e) list())
    curInnerOpt <- tryCatch(baseCtl$innerOpt, error = function(e) NA_integer_)
    if (identical(curInnerOpt, 1L)) {
      n1 <- list(ok = TRUE, time = unname(t1 - t0),
                 objf = tryCatch(as.numeric(fit$objf), error = function(e) NA_real_),
                 theta = tryCatch(setNames(as.numeric(fit$theta), names(fit$theta)), error = function(e) NULL),
                 error = NA_character_)
    } else {
      n1 <- .benchFitOnce(object, data, est, baseCtl, table, dots, save, envir, "n1qn1")
    }
    tr <- .benchFitOnce(object, data, est, baseCtl, table, dots, save, envir, "trust")
    .benchRecord(.benchSourceTag, label, n1, tr)
  }
  fit
}

.withBenchOverrides <- function(envir) {
  assign("nlmixr2", .benchNlmixr2Wrap, envir = envir)
  assign("nlmixr", .benchNlmixr2Wrap, envir = envir)
  # Neutralize testthat's skip/report machinery -- run the block, don't assert.
  assign("test_that", function(desc, code) {
    tryCatch(force(code), error = function(e) NULL)
  }, envir = envir)
  assign("nmTest", function(code) {
    tryCatch(force(code), error = function(e) NULL)
  }, envir = envir)
  assign("skip_on_cran", function(...) invisible(NULL), envir = envir)
  assign("skip_if_not_installed", function(...) invisible(NULL), envir = envir)
  assign("skip_on_ci", function(...) invisible(NULL), envir = envir)
  assign("skip", function(...) invisible(NULL), envir = envir)
}

.sourceCorpusFile <- function(path, sourceTag) {
  .benchSourceTag <<- sourceTag
  .benchModelCounter <<- 0L
  env <- new.env(parent = globalenv())
  .withBenchOverrides(env)
  ok <- tryCatch({
    source(path, local = env, echo = FALSE)
    TRUE
  }, error = function(e) {
    .blog("[%s] source() failed: %s", sourceTag, conditionMessage(e))
    FALSE
  })
  if (ok && .benchModelCounter == 0L) {
    .benchLogSkip(sourceTag, "(whole file)", "sourced cleanly but no focei fit was intercepted")
  } else if (!ok) {
    .benchLogSkip(sourceTag, "(whole file)", "source() failed -- see log")
  }
  invisible(ok)
}

# ---------------------------------------------------------------------------
# (a) discussion #924's 5 benchmark models, hand-written directly (no
# interception needed -- run explicitly at each innerOpt).
# ---------------------------------------------------------------------------
.runDiscussion924 <- function() {
  .theoOneCmt <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  .theoBlockOmega <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6
      eta.cl + eta.v ~ c(0.3, 0.01, 0.1)
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  .oral2cmt <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; tv2 <- 3; tq <- 0.5
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      v2 <- exp(tv2); q <- exp(tq)
      linCmt() ~ add(add.sd)
    })
  }
  .oral1cmtMM <- function() {
    ini({
      tka <- 0.45; tvm <- 2; tkm <- 3; tv <- 3.45
      eta.ka ~ 0.6; eta.vm ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); vm <- exp(tvm + eta.vm); km <- exp(tkm); v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - (vm / (km + center / v)) * (center / v)
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  .phenoSparse <- function() {
    ini({
      tcl <- 0.02; tv <- 1.5
      eta.cl ~ 0.1; eta.v ~ 0.1
      add.sd <- 0.5
    })
    model({
      cl <- tcl * exp(eta.cl); v <- tv * exp(eta.v)
      cp <- center / v
      d/dt(center) <- -cl / v * center
      cp ~ add(add.sd)
    })
  }
  cases <- list(
    list(label = "d924-theo-1cmt", ui = .theoOneCmt, data = nlmixr2data::theo_sd),
    list(label = "d924-theo-block-omega", ui = .theoBlockOmega, data = nlmixr2data::theo_sd),
    list(label = "d924-oral-2cmt", ui = .oral2cmt, data = nlmixr2data::theo_sd),
    list(label = "d924-oral-1cmt-mm", ui = .oral1cmtMM, data = nlmixr2data::theo_sd),
    list(label = "d924-pheno-sparse", ui = .phenoSparse, data = nlmixr2data::pheno_sd)
  )
  for (cs in cases) {
    ok <- tryCatch({
      n1 <- .benchFitOnce(cs$ui, cs$data, "focei", list(), nlmixr2est::tableControl(), list(),
                           NULL, parent.frame(), "n1qn1")
      tr <- .benchFitOnce(cs$ui, cs$data, "focei", list(), nlmixr2est::tableControl(), list(),
                           NULL, parent.frame(), "trust")
      .benchRecord("discussion-924", cs$label, n1, tr)
      TRUE
    }, error = function(e) { .blog("[discussion-924] %s errored: %s", cs$label, conditionMessage(e)); FALSE })
    if (!ok) .benchLogSkip("discussion-924", cs$label, "case setup errored -- see log")
  }
}

# ---------------------------------------------------------------------------
# (b) nlmixr2 package vignettes
# ---------------------------------------------------------------------------
.runVignettes <- function(vigDir) {
  if (!dir.exists(vigDir)) {
    .blog("vignette directory not found: %s -- skipping corpus (b)", vigDir)
    return(invisible(NULL))
  }
  if (!requireNamespace("knitr", quietly = TRUE)) {
    .blog("knitr not installed -- skipping vignette corpus (b)")
    return(invisible(NULL))
  }
  rmds <- sort(list.files(vigDir, pattern = "\\.Rmd$", recursive = TRUE, full.names = TRUE))
  .blog("found %d vignette .Rmd files under %s", length(rmds), vigDir)
  for (rmd in rmds) {
    tag <- paste0("vignette:", sub("\\.Rmd$", "", basename(rmd)))
    rFile <- tempfile(fileext = ".R")
    purlOk <- tryCatch({
      knitr::purl(rmd, output = rFile, documentation = 0L, quiet = TRUE)
      TRUE
    }, error = function(e) {
      .blog("[%s] knitr::purl() failed: %s", tag, conditionMessage(e))
      FALSE
    })
    if (!purlOk) { .benchLogSkip(tag, "(whole file)", "knitr::purl() failed -- see log"); next }
    .sourceCorpusFile(rFile, tag)
    unlink(rFile)
  }
}

# ---------------------------------------------------------------------------
# (c) non-wang2007 tests/testthat/test-focei-*.R
# ---------------------------------------------------------------------------
.runFoceiTests <- function(testDir = "tests/testthat") {
  if (!dir.exists(testDir)) {
    .blog("test directory not found: %s -- skipping corpus (c)", testDir)
    return(invisible(NULL))
  }
  files <- sort(list.files(testDir, pattern = "^test-focei-.*\\.R$", full.names = TRUE))
  files <- files[!grepl("wang2007", basename(files))]
  # This harness's own test file benchmarks nothing new -- exclude it.
  files <- files[basename(files) != "test-focei-trust-inner.R"]
  .blog("found %d non-wang2007 test-focei-*.R files under %s", length(files), testDir)
  for (f in files) {
    tag <- paste0("test:", sub("\\.R$", "", basename(f)))
    .sourceCorpusFile(f, tag)
  }
}

# ---------------------------------------------------------------------------
# Run all three corpora, then write a markdown summary.
# ---------------------------------------------------------------------------
.blog("--- (a) discussion #924 models ---")
.runDiscussion924()
.blog("--- (b) nlmixr2 vignettes (%s) ---", .vigDir)
.runVignettes(.vigDir)
.blog("--- (c) non-wang2007 test-focei-*.R ---")
.runFoceiTests()

if (length(.benchResults)) {
  df <- do.call(rbind, .benchResults)
  write.csv(df, .csvPath, row.names = FALSE)
  mdPath <- file.path(.outDir, "trust-inner-benchmark.md")
  ok <- df[!is.na(df$n1qn1_ok) & df$n1qn1_ok & !is.na(df$trust_ok) & df$trust_ok, , drop = FALSE]
  lines <- c(
    "# innerOpt=\"trust\" vs innerOpt=\"n1qn1\" benchmark",
    "",
    sprintf("Total corpus entries: %d (both converged: %d, skipped/failed: %d)",
            nrow(df), nrow(ok), nrow(df) - nrow(ok)),
    ""
  )
  if (nrow(ok)) {
    lines <- c(lines,
      sprintf("Median speedup (n1qn1 time / trust time): %.2fx", stats::median(ok$speedup, na.rm = TRUE)),
      sprintf("Median |objf diff|: %.4g", stats::median(abs(ok$dObjf), na.rm = TRUE)),
      sprintf("Median max |param diff|: %.4g", stats::median(ok$maxParDiff, na.rm = TRUE)),
      ""
    )
  }
  lines <- c(lines, "| source | model | n1qn1 ok | n1qn1 time | n1qn1 objf | trust ok | trust time | trust objf | dObjf | speedup | maxParDiff |",
             "|---|---|---|---|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    lines <- c(lines, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
      r$source, r$model, r$n1qn1_ok, format(r$n1qn1_time), format(r$n1qn1_objf),
      r$trust_ok, format(r$trust_time), format(r$trust_objf),
      format(r$dObjf), format(r$speedup), format(r$maxParDiff)))
  }
  writeLines(lines, mdPath)
  .blog("wrote %s and %s", .csvPath, mdPath)
} else {
  .blog("no results recorded at all -- something is wrong with the harness or the corpora")
}
close(.logCon)
