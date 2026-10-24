#!/usr/bin/env Rscript
# Corpus benchmark of the FOCEi-family inner (eta) optimizers, the #1153 corpus:
# every test file that fits a FOCEi-family method, the nlmixr2 vignettes and the
# discussion #924 models, each case from 3 deterministic starts.
#
# Usage (run from the package root; pkgload loads this checkout):
#   Rscript inst/benchmarks/benchmark-inner-opt.R capture <outDir> [vignetteDir]
#   Rscript inst/benchmarks/benchmark-inner-opt.R replay  <outDir> <arm> [threads]
#   Rscript inst/benchmarks/benchmark-inner-opt.R report  <outDir> [refArm]
#
# capture sources the corpus with nlmixr2()/nlmixr() intercepted and saves each
# FOCEi-family call (model, data, est, control) to <outDir>/cases/.  replay
# re-runs every case under one arm (see .arms) and appends to
# <outDir>/arm-<arm>.csv; arms can run concurrently.  report compares the arms.

.args <- commandArgs(trailingOnly = TRUE)
.mode <- .args[1]
.outDir <- .args[2]
dir.create(file.path(.outDir, "cases"), showWarnings = FALSE, recursive = TRUE)

# FOCEi-family estimation methods (they run the inner eta problem)
.foceiFamily <- c("focei", "foce", "foi", "posthoc", "agq", "laplace", "focep",
                  "mfocei", "mfoce", "mfocep", "mlaplace", "magq", "impmap")

# arm -> control fields set on every case
.arms <- list(
  trust = list(innerOpt = 3L, trustPolish = FALSE),
  trustPolish = list(innerOpt = 3L, trustPolish = TRUE),
  steihaug = list(innerOpt = 6L, steihaugPolish = FALSE),
  steihaugPolish = list(innerOpt = 6L, steihaugPolish = TRUE),
  newuoa = list(innerOpt = 5L, newuoaPolish = FALSE),
  newuoaPolish = list(innerOpt = 5L, newuoaPolish = TRUE)
)

.loadPkg <- function() {
  suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE, helpers = FALSE))
}

# ---------------------------------------------------------------------------
# capture
# ---------------------------------------------------------------------------
.capture <- function(vigDir) {
  .loadPkg()
  rxode2::setRxThreads(4L)
  .logCon <- file(file.path(.outDir, "capture.log"), open = "a")
  .clog <- function(fmt, ...) {
    cat(format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...), "\n", file = .logCon)
    flush(.logCon)
  }
  .n <- 0L
  .tag <- ""
  .orig <- nlmixr2est::nlmixr2
  .wrap <- function(object, data, est = NULL, control = list(),
                    table = nlmixr2est::tableControl(), ..., save = NULL,
                    envir = parent.frame()) {
    .est <- if (is.null(est)) NA_character_ else as.character(est)
    if (!is.na(.est) && .est %in% .foceiFamily && !missing(data)) {
      .n <<- .n + 1L
      .case <- try(list(source = .tag, idx = .n, est = .est,
                        ui = rxode2::rxUiDecompress(rxode2::assertRxUi(object)),
                        data = as.data.frame(data), control = control,
                        dots = list(...)), silent = TRUE)
      if (inherits(.case, "try-error")) {
        .clog("[%s] #%d not captured: %s", .tag, .n, .case)
      } else {
        saveRDS(.case, file.path(.outDir, "cases",
                                 sprintf("%s__%03d.rds", gsub("[^A-Za-z0-9_-]", "_", .tag), .n)))
      }
    }
    .orig(object = object, data = data, est = est, control = control, table = table,
          ..., save = save, envir = envir)
  }
  .base <- new.env(parent = asNamespace("nlmixr2est"))
  assign("nlmixr2", .wrap, envir = .base)
  assign("nlmixr", .wrap, envir = .base)
  assign("test_that", function(desc, code) tryCatch(force(code), error = function(e) NULL), envir = .base)
  assign("nmTest", function(code, ...) tryCatch(force(code), error = function(e) NULL), envir = .base)
  for (.s in c("skip_on_cran", "skip_if_not_installed", "skip_on_ci", "skip",
               "skip_if", "skip_if_not", "skip_on_os")) {
    assign(.s, function(...) invisible(NULL), envir = .base)
  }
  suppressPackageStartupMessages(library(testthat))
  for (.h in list.files("tests/testthat", "^helper-.*\\.R$", full.names = TRUE)) {
    if (grepl("zzz-fits", .h)) next
    try(sys.source(.h, envir = .base), silent = TRUE)
  }
  .runFile <- function(path, tag) {
    .tag <<- tag
    .n <<- 0L
    .env <- new.env(parent = .base)
    .t0 <- proc.time()[["elapsed"]]
    .ok <- try(suppressWarnings(suppressMessages(sys.source(path, envir = .env))), silent = TRUE)
    .clog("[%s] %d captured in %.0fs%s", tag, .n, proc.time()[["elapsed"]] - .t0,
          if (inherits(.ok, "try-error")) paste0(" (source error: ", .ok, ")") else "")
  }
  # test files that fit a FOCEi-family method
  .re <- paste0("est *= *\"(", paste(.foceiFamily, collapse = "|"), ")\"")
  for (.f in sort(list.files("tests/testthat", "^test-.*\\.R$", full.names = TRUE))) {
    if (!any(grepl(.re, readLines(.f, warn = FALSE)))) next
    .runFile(.f, paste0("test_", sub("^test-(.*)\\.R$", "\\1", basename(.f))))
  }
  # vignettes
  for (.rmd in sort(list.files(vigDir, "\\.Rmd$", full.names = TRUE))) {
    .r <- tempfile(fileext = ".R")
    if (inherits(try(knitr::purl(.rmd, output = .r, documentation = 0L, quiet = TRUE),
                     silent = TRUE), "try-error")) next
    .old <- setwd(dirname(.rmd))
    .runFile(.r, paste0("vig_", sub("\\.Rmd$", "", basename(.rmd))))
    setwd(.old)
  }
  # discussion #924 models
  .d924 <- file.path(.outDir, "d924.R")
  writeLines(c(
    "f <- function(m, d) nlmixr(m, d, est = 'focei', control = foceiControl(print = 0))",
    "f(function() { ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })",
    "  model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v); linCmt() ~ add(add.sd) }) }, nlmixr2data::theo_sd)",
    "f(function() { ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; eta.cl + eta.v ~ c(0.3, 0.01, 0.1); add.sd <- 0.7 })",
    "  model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v); linCmt() ~ add(add.sd) }) }, nlmixr2data::theo_sd)",
    "f(function() { ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; tv2 <- 3; tq <- 0.5; eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })",
    "  model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v); v2 <- exp(tv2); q <- exp(tq); linCmt() ~ add(add.sd) }) }, nlmixr2data::theo_sd)",
    "f(function() { ini({ tka <- 0.45; tvm <- 2; tkm <- 3; tv <- 3.45; eta.ka ~ 0.6; eta.vm ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })",
    "  model({ ka <- exp(tka + eta.ka); vm <- exp(tvm + eta.vm); km <- exp(tkm); v <- exp(tv + eta.v)",
    "    d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - (vm / (km + center / v)) * (center / v)",
    "    cp <- center / v; cp ~ add(add.sd) }) }, nlmixr2data::theo_sd)",
    "f(function() { ini({ tcl <- 0.02; tv <- 1.5; eta.cl ~ 0.1; eta.v ~ 0.1; add.sd <- 0.5 })",
    "  model({ cl <- tcl * exp(eta.cl); v <- tv * exp(eta.v); cp <- center / v; d/dt(center) <- -cl / v * center; cp ~ add(add.sd) }) }, nlmixr2data::pheno_sd)"
  ), .d924)
  .runFile(.d924, "d924")
  close(.logCon)
}

# ---------------------------------------------------------------------------
# replay
# ---------------------------------------------------------------------------
# start 0: as captured; 1: thetas +10%, omega x1.25; 2: thetas -10%, omega x0.8
.startUi <- function(ui, start) {
  if (start == 0L) return(ui)
  .df <- ui$iniDf
  .th <- !is.na(.df$ntheta) & !.df$fix
  .f <- if (start == 1L) 1.1 else 0.9
  .df$est[.th] <- pmin(pmax(.df$est[.th] * .f, .df$lower[.th]), .df$upper[.th])
  .om <- !is.na(.df$neta1) & !.df$fix
  .df$est[.om] <- .df$est[.om] * (if (start == 1L) 1.25 else 0.8)
  rxode2::ini(ui) <- .df
  ui
}

.features <- function(ui, data, control, est) {
  .pd <- ui$predDf
  .lst <- paste(deparse(ui$lstExpr), collapse = " ")
  .df <- ui$iniDf
  .om <- .df[!is.na(.df$neta1), ]
  list(
    errAdd = any(.pd$errType == "add" & .pd$errTypeF == "none"),
    errProp = any(.pd$errType == "prop"),
    errAddProp = any(.pd$errType == "add + prop"),
    linCmt = grepl("linCmt", .lst),
    ode = grepl("d/dt", .lst, fixed = TRUE),
    blockOmega = any(.om$neta1 != .om$neta2),
    multiEndpoint = nrow(.pd) > 1L,
    cens = any(toupper(names(data)) == "CENS") && any(data[[which(toupper(names(data)) == "CENS")[1]]] != 0, na.rm = TRUE),
    covariates = length(ui$allCovs) > 0L,
    fast = isTRUE(control$fast),
    eta3 = sum(.om$neta1 == .om$neta2) >= 3L,
    est = est
  )
}

.noFeatures <- function(est) {
  list(errAdd = NA, errProp = NA, errAddProp = NA, linCmt = NA, ode = NA,
       blockOmega = NA, multiEndpoint = NA, cens = NA, covariates = NA,
       fast = NA, eta3 = NA, est = est)
}

.replay <- function(arm, threads) {
  .loadPkg()
  rxode2::setRxThreads(threads)
  .set <- .arms[[arm]]
  if (is.null(.set)) stop("unknown arm ", arm)
  .csv <- file.path(.outDir, paste0("arm-", arm, ".csv"))
  .done <- if (file.exists(.csv)) unique(paste(utils::read.csv(.csv)$case, utils::read.csv(.csv)$start)) else character(0)
  for (.p in sort(list.files(file.path(.outDir, "cases"), "\\.rds$", full.names = TRUE))) {
    .c <- readRDS(.p)
    .case <- sub("\\.rds$", "", basename(.p))
    for (.start in 0:2) {
      if (paste(.case, .start) %in% .done) next
      .ctl <- .c$control
      if (is.null(.ctl)) .ctl <- list()
      for (.n in names(.set)) .ctl[[.n]] <- .set[[.n]]
      if (is.list(.ctl) && !inherits(.ctl, "foceiControl") && is.null(attr(.ctl, "class"))) {
        # an unbuilt list: hand the method its names
        .iof <- utils::getFromNamespace(".innerOptFun", "nlmixr2est")
        .ctl$innerOpt <- names(.iof)[match(.set$innerOpt, .iof)]
      }
      .ctl$print <- 0L
      .ctl$calcTables <- FALSE
      .row <- list(case = .case, source = .c$source, start = .start, arm = arm,
                   ok = FALSE, objf = NA_real_, time = NA_real_, maxOuter = NA_integer_,
                   calls = NA_integer_, polish = NA_integer_, failed = NA_integer_,
                   nRunInfo = NA_integer_, error = NA_character_)
      .ui <- try(.startUi(.c$ui, .start), silent = TRUE)
      .t0 <- proc.time()[["elapsed"]]
      .fit <- if (inherits(.ui, "try-error")) .ui else
        try(suppressWarnings(suppressMessages(do.call(nlmixr2est::nlmixr2,
          c(list(object = .ui, data = .c$data, est = .c$est, control = .ctl), .c$dots)))),
          silent = TRUE)
      .row$time <- proc.time()[["elapsed"]] - .t0
      if (inherits(.fit, "try-error")) {
        .row$error <- substr(gsub("[\r\n]+", " ", as.character(.fit)), 1, 200)
      } else {
        .row$ok <- TRUE
        .row$objf <- tryCatch(as.numeric(.fit$objf), error = function(e) NA_real_)
        .row$maxOuter <- tryCatch(as.integer(.fit$control$maxOuterIterations), error = function(e) NA_integer_)
        .nt <- tryCatch(.fit$env$nTrustInner, error = function(e) NULL)
        if (!is.null(.nt)) {
          .row$calls <- .nt[["calls"]]
          if ("polish" %in% names(.nt)) .row$polish <- .nt[["polish"]]
          if ("failed" %in% names(.nt)) .row$failed <- .nt[["failed"]]
        }
        .row$nRunInfo <- tryCatch(length(.fit$runInfo), error = function(e) NA_integer_)
      }
      .ft <- try(.features(.c$ui, .c$data, .c$control, .c$est), silent = TRUE)
      if (inherits(.ft, "try-error")) .ft <- .noFeatures(.c$est)
      .row <- c(.row, .ft)
      .df <- as.data.frame(.row, stringsAsFactors = FALSE)
      utils::write.table(.df, .csv, sep = ",", row.names = FALSE,
                         col.names = !file.exists(.csv), append = file.exists(.csv))
    }
  }
}

if (identical(.mode, "capture")) {
  .capture(if (length(.args) >= 3) .args[3] else path.expand("~/src/nlmixr2/vignettes"))
} else if (identical(.mode, "replay")) {
  .replay(.args[3], if (length(.args) >= 4) as.integer(.args[4]) else 4L)
} else if (identical(.mode, "report")) {
  source(file.path("inst", "benchmarks", "benchmark-inner-opt-report.R"))
  .innerOptReport(.outDir, if (length(.args) >= 3) .args[3] else "trust")
} else {
  stop("mode must be capture, replay or report")
}
