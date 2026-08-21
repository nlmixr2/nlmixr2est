#!/usr/bin/env Rscript
# Benchmark foceiControl(hessianMethod=) -- the quasi-Newton update for
# calcEtaHessian()'s non-normal-endpoint (needOptimHess) inner Hessian
# (src/inner.cpp) -- against the default "fd" (finite-difference-of-the-
# gradient every Newton step) on a handful of representative non-normal-
# endpoint models. This is a DIFFERENT code path from
# inst/benchmarks/benchmark-trust-outer.R's est="trust" outer-theta sweep:
# normal-endpoint models never exercise this branch at all (the Gauss-Newton
# inner Hessian is used unconditionally there), so this benchmark is built
# from small, purpose-written non-normal-endpoint models rather than this
# package's existing (mostly normal-endpoint) nlm-family test corpus.
#
# Usage:
#   Rscript inst/benchmarks/benchmark-focei-hessian-method.R [outDir]
# outDir defaults to inst/benchmarks/results.

suppressPackageStartupMessages({
  library(nlmixr2est)
})

.bargs <- commandArgs(trailingOnly = TRUE)
.outDir <- if (length(.bargs) >= 1) .bargs[1] else "inst/benchmarks/results"
dir.create(.outDir, showWarnings = FALSE, recursive = TRUE)
.csvPath <- file.path(.outDir, "focei-hessian-method-benchmark.csv")
.logPath <- file.path(.outDir, "focei-hessian-method-benchmark.log")
.mdPath  <- file.path(.outDir, "focei-hessian-method-benchmark.md")

.logCon <- file(.logPath, open = "a")
.blog <- function(fmt, ...) {
  msg <- sprintf(fmt, ...)
  cat(msg, "\n")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg, "\n", file = .logCon)
  flush(.logCon)
}
.blog("=== benchmark-focei-hessian-method.R starting ===")

.hessianMethods <- c("fd", "bfgs", "sr1", "bofill")

.benchResults <- list()
.record <- function(model, results) {
  row <- data.frame(model = model, stringsAsFactors = FALSE)
  fdObjf <- results[["fd"]]$objf
  for (hm in .hessianMethods) {
    r <- results[[hm]]
    row[[paste0(hm, "_ok")]] <- isTRUE(r$ok)
    row[[paste0(hm, "_time")]] <- as.numeric(r$time)
    row[[paste0(hm, "_objf")]] <- as.numeric(r$objf)
    row[[paste0(hm, "_dObjf")]] <- as.numeric(r$objf) - as.numeric(fdObjf)
    row[[paste0(hm, "_nHessianQN")]] <- as.integer(r$nHessianQN)
    row[[paste0(hm, "_error")]] <- ifelse(is.null(r$error), NA_character_, r$error)
  }
  .benchResults[[length(.benchResults) + 1]] <<- row
  write.csv(do.call(rbind, .benchResults), .csvPath, row.names = FALSE)
  .blog("[%s] %s", model,
        paste(sapply(.hessianMethods, function(hm) {
          r <- results[[hm]]
          sprintf("%s: ok=%s time=%.2fs objf=%s nHessianQN=%d", hm, r$ok, r$time,
                  format(r$objf), r$nHessianQN)
        }), collapse = " | "))
}

.fitOnce <- function(object, data, hessianMethod, ...) {
  t0 <- proc.time()["elapsed"]
  fit <- tryCatch(
    suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(object, data, est = "focei",
                          control = nlmixr2est::foceiControl(print = 0L, hessianMethod = hessianMethod, ...))
    )),
    error = function(e) e
  )
  t1 <- proc.time()["elapsed"]
  if (inherits(fit, "error")) {
    return(list(ok = FALSE, time = unname(t1 - t0), objf = NA_real_,
                nHessianQN = nlmixr2est:::.nHessianQN(), error = conditionMessage(fit)))
  }
  list(ok = TRUE, time = unname(t1 - t0),
       objf = tryCatch(as.numeric(fit$objective), error = function(e) NA_real_),
       nHessianQN = nlmixr2est:::.nHessianQN(), error = NA_character_)
}

.runModel <- function(label, object, data, ...) {
  results <- list()
  for (hm in .hessianMethods) {
    results[[hm]] <- .fitOnce(object, data, hm, ...)
  }
  .record(label, results)
}

# ---------------------------------------------------------------------------
# Corpus: small, fast-converging non-normal-endpoint models.
# ---------------------------------------------------------------------------

set.seed(1)

# 1. Poisson count model.
poisMod <- function() {
  ini({
    te0 <- log(5)
    eta.e0 ~ 0.3
  })
  model({
    e0 <- exp(te0 + eta.e0)
    effect <- e0
    effect ~ dpois(effect)
  })
}
poisData <- data.frame(ID = rep(1:30, each = 4), TIME = rep(c(0, 1, 2, 3), 30))
poisData$DV <- rpois(nrow(poisData), lambda = 5 * exp(rnorm(30, 0, 0.3))[poisData$ID])
.runModel("poisson-e0", poisMod, poisData)

# NOTE: a binomial (dbinom) response model was tried here and dropped -- it
# crashes the R session (Rcpp::LongjumpException) even under the DEFAULT
# hessianMethod="fd", i.e. with zero of this feature's code engaged. That is
# a pre-existing issue elsewhere (or in how that toy model was specified),
# not something hessianMethod= causes; not chased further here (out of scope).

# 2. General ll() Gaussian-equivalent (sanity check: should match an
# equivalent add()-error fit's objective closely, since ll() here is just
# the Gaussian log-density spelled out explicitly).
llMod <- function() {
  ini({
    tka <- 0.45; tcl <- 1; tv <- 3.45
    eta.ka ~ 0.3; eta.cl ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
    d/dt(depot) <- -ka * depot
    d/dt(centr) <- ka * depot - cl / v * centr
    cp <- centr / v
    ll(err) ~ -log(add.sd) - 0.5*log(2*pi) - 0.5*((DV-cp)/add.sd)^2
  })
}
llData <- nlmixr2data::theo_sd
.runModel("ll-gaussian-equivalent", llMod, llData)

# ---------------------------------------------------------------------------
# Markdown summary.
# ---------------------------------------------------------------------------
if (length(.benchResults)) {
  df <- do.call(rbind, .benchResults)
  write.csv(df, .csvPath, row.names = FALSE)
  lines <- c(
    "# foceiControl(hessianMethod=) benchmark (non-normal-endpoint models)",
    "",
    sprintf("Total models: %d", nrow(df)),
    ""
  )
  for (hm in .hessianMethods) {
    okCol <- paste0(hm, "_ok")
    ok <- df[!is.na(df[[okCol]]) & df[[okCol]], , drop = FALSE]
    lines <- c(lines, sprintf("## hessianMethod=\"%s\" (converged: %d/%d)", hm, nrow(ok), nrow(df)), "")
    if (nrow(ok)) {
      lines <- c(lines,
        sprintf("- Median time: %.3fs", stats::median(ok[[paste0(hm, "_time")]], na.rm = TRUE)),
        sprintf("- Median |objf diff| vs fd: %.4g", stats::median(abs(ok[[paste0(hm, "_dObjf")]]), na.rm = TRUE)),
        sprintf("- Total nHessianQN calls: %d", sum(ok[[paste0(hm, "_nHessianQN")]], na.rm = TRUE)),
        "")
    }
  }
  header <- paste0("| model | ", paste(sapply(.hessianMethods, function(hm)
    sprintf("%s ok | %s time | %s objf | %s dObjf", hm, hm, hm, hm)), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", 1 + 4 * length(.hessianMethods)), collapse = "|"), "|")
  lines <- c(lines, header, sep)
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    cells <- r$model
    for (hm in .hessianMethods) {
      cells <- c(cells, r[[paste0(hm, "_ok")]], format(r[[paste0(hm, "_time")]]),
                 format(r[[paste0(hm, "_objf")]]), format(r[[paste0(hm, "_dObjf")]]))
    }
    lines <- c(lines, paste0("| ", paste(cells, collapse = " | "), " |"))
  }
  writeLines(lines, .mdPath)
  .blog("wrote %s and %s", .csvPath, .mdPath)
} else {
  .blog("no results recorded -- something is wrong with the harness or corpus")
}
close(.logCon)
