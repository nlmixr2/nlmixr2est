## Phase 4.1: the four Bauer arms against the four estimators, one pass.
##
## The arms are Bauer's own datasets and control files (gamma_indpar/), and they
## differ ONLY in the declared relative variance -- theta 5 and 6 of each
## _sim.ctl.  That is the whole point of the set: the same PK model and the same
## design at four dispersions, so a method that is fine at CV 30% and falls over
## at CV 141% shows it here and nowhere else.
##
##   g1  lclrv -2.408     rv 0.09   CV  30%
##   g3  lclrv -0.693147  rv 0.5    CV  71%
##   g2  lclrv  0.00001   rv 1.0    CV 100%
##   g4  lclrv  0.693147  rv 2.0    CV 141%
##
## Starting values are Bauer's own (_imp.ctl $THETA), displaced from the truth
## on every arm, so "the fit recovers it" cannot pass by inertia.
##
## Run:  Rscript benchArms.R <arm> <method>   (one cell, so a cell that hangs
## costs one cell)  or  Rscript benchArms.R table   to print what has finished.
##
## `mceta` is held FIXED and reported in every row: imp inherits its effect
## through the shared FOCEi inner MAP, and it is recorded as blowing up one
## iteration of exactly this model (src/imp.cpp), so a table that let it vary
## would not be comparing methods.
.a <- commandArgs(TRUE)
D <- Sys.getenv("BENCH_DIR", file.path(tempdir(), "benchArms"))
dir.create(D, showWarnings = FALSE, recursive = TRUE)
.lib <- Sys.getenv("BENCH_LIB", "")
if (nzchar(.lib)) .libPaths(c(.lib, .libPaths()))

DAT <- Sys.getenv("BENCH_DAT", "~/src/gamma_indpar")
MCETA <- as.integer(Sys.getenv("BENCH_MCETA", "0"))

## truth (from <arm>_sim.ctl) and starts (from <arm>_imp.ctl), transcribed so
## the table is readable without the NONMEM files to hand
ARMS <- list(
  g1 = list(file = "gamma_clv1.dat",  rv = -2.408,    start = -3.0,  rho0 = 0.3),
  g3 = list(file = "gamma3_clv1.dat", rv = -0.693147, start = -0.6,  rho0 = 0.3),
  g2 = list(file = "gamma2_clv1.dat", rv = 0.00001,   start = -0.1,  rho0 = 0.3),
  g4 = list(file = "gamma4_clv1.dat", rv = 0.693147,  start = 0.7,   rho0 = 0.6)
)
TRUTH <- c(lclm = 1.63, lv1m = 1.55, tq = 0.74, tv2 = 2.3, rho = 0.5)

.mkModel <- function(arm) {
  .s <- ARMS[[arm]]
  eval(parse(text = sprintf('function() {
    ini({
      lclm <- 1.9; lv1m <- 1.8; tq <- 0.9; tv2 <- 4.2
      lclrv <- %s; lv1rv <- %s
      eta.cl + eta.v1 ~ c(1, %s, 1)
      dist(eta.cl) ~ dgamma(shape = 1 / exp(lclrv),
                            rate = 1 / (exp(lclrv) * exp(lclm)))
      dist(eta.v1) ~ dgamma(shape = 1 / exp(lv1rv),
                            rate = 1 / (exp(lv1rv) * exp(lv1m)))
      eta.q + eta.v2 ~ c(0.1, 0.01, 0.1)
      prop.sd <- 0.316
    })
    model({
      cl <- eta.cl; v <- eta.v1
      q <- exp(tq + eta.q); v2 <- exp(tv2 + eta.v2)
      linCmt() ~ prop(prop.sd)
    })
  }', .s$start, .s$start, .s$rho0)))
}

.control <- function(method) {
  switch(method,
    saem  = nlmixr2est::saemControl(print = 0L, covMethod = "", nBurn = 300L, nEm = 300L),
    focei = nlmixr2est::foceiControl(print = 0L, covMethod = "", mceta = MCETA),
    imp   = nlmixr2est::impmapControl(print = 0L, covMethod = "", nIter = 100L,
                                      mapIter = 0L, mceta = MCETA),
    vae   = nlmixr2est::vaeControl(print = 0L, covMethod = ""),
    stop("unknown method: ", method))
}

## Two MAREs, not one.  A single mean over means AND relative variances is not
## comparable ACROSS arms, which is the one comparison this table exists to
## make: the arms differ only in rv, so the rv denominator moves by 22x from g1
## (0.09) to g4 (2.0).  The same absolute miss on rv is ~490% on g1 and ~20% on
## g4, and a combined MARE would report g1 as catastrophic and g4 as fine for
## identical accuracy.  Reported separately, each column is comparable down its
## own arm and across arms.
##
## Compared as what they MEAN -- exp() of the log-mean is a clearance, exp() of
## the log relative variance is a relative variance -- because a relative error
## on a log scale is not a relative error on the quantity.
.mareParts <- function(p, arm) {
  .tm <- c(exp(TRUTH[["lclm"]]), exp(TRUTH[["lv1m"]]))
  .em <- c(exp(p[["lclm"]]), exp(p[["lv1m"]]))
  .tv <- rep(exp(ARMS[[arm]]$rv), 2)
  .ev <- c(exp(p[["lclrv"]]), exp(p[["lv1rv"]]))
  c(mean = 100 * mean(abs(.em - .tm) / .tm),
    rv   = 100 * mean(abs(.ev - .tv) / .tv))
}

if (identical(.a[1], "table")) {
  fs <- list.files(D, "^cell-.*[.]rds$", full.names = TRUE)
  if (!length(fs)) { cat("no cells finished yet\n"); quit(status = 0) }
  rows <- lapply(fs, function(f) {
    r <- readRDS(f)
    data.frame(arm = r$arm, cv = sprintf("%.0f%%", 100 * sqrt(exp(ARMS[[r$arm]]$rv))),
               method = r$method,
               mareMean = if (is.null(r$mare)) NA_real_ else round(r$mare[["mean"]], 1),
               mareRv = if (is.null(r$mare)) NA_real_ else round(r$mare[["rv"]], 1),
               rho = if (is.null(r$p) || !("rxCor.eta.v1.eta.cl" %in% names(r$p))) NA_real_
                     else round(tanh(r$p[["rxCor.eta.v1.eta.cl"]]), 3),
               objf = if (is.null(r$objf)) NA_real_ else round(r$objf, 1),
               secs = round(r$el), mceta = r$mceta,
               note = if (is.null(r$err)) "" else substr(r$err, 1, 40),
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(match(tab$arm, c("g1","g3","g2","g4")), tab$method), ]
  print(tab, row.names = FALSE)
  quit(status = 0)
}

arm <- .a[1]; method <- .a[2]
stopifnot(arm %in% names(ARMS))
suppressMessages(library(nlmixr2est)); suppressMessages(library(rxode2))
rxode2::setRxThreads(as.integer(Sys.getenv("BENCH_CORES", "5")))
d <- read.table(file.path(path.expand(DAT), ARMS[[arm]]$file), skip = 1, header = TRUE)
names(d) <- toupper(names(d))
t0 <- proc.time()[["elapsed"]]
f <- try(suppressWarnings(nlmixr2(.mkModel(arm), d, est = method,
                                  control = .control(method))), silent = TRUE)
el <- proc.time()[["elapsed"]] - t0
out <- list(arm = arm, method = method, el = el, mceta = MCETA)
if (inherits(f, "try-error")) {
  out$err <- conditionMessage(attr(f, "condition"))
  cat("CELL", arm, method, "ERROR after", round(el), "s\n")
} else {
  p <- setNames(f$parFixedDf$Estimate, rownames(f$parFixedDf))
  out$p <- p; out$objf <- as.numeric(f$objf); out$mare <- .mareParts(p, arm)
  cat("CELL", arm, method, "OK", round(el), "s  MARE mean",
      round(out$mare[["mean"]], 1), "% rv", round(out$mare[["rv"]], 1),
      "%  objf", round(out$objf, 1), "\n")
}
saveRDS(out, file.path(D, sprintf("cell-%s-%s.rds", arm, method)))
