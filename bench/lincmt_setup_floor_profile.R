# lincmt_setup_floor_profile.R -- where does the ~1.2-1.4 s per-fit setup
# floor of a linCmt() FOCEi posthoc fit go?  (Companion to
# rxode2-lincmt-carry-jump/bench/lincmt_vs_ode_focei.R, which measured the
# floor; this script attributes it stage by stage.)
#
# Protocol: optimized rxode2 (~/src/rxode2-lincmt-carry-jump, compile=FALSE),
# this nlmixr2est worktree (helpers=FALSE), single rxode2 thread, run pinned:
#   REPS=3 taskset -c <idle core> Rscript bench/lincmt_setup_floor_profile.R
# Reports: (1) wall time of rep 1 (cold: fresh session, model text never seen)
# vs warm repeats of the identical fit; (2) an Rprof by.total breakdown of one
# warm repeat, aggregated to stages by function name.

suppressMessages({
  devtools::load_all("~/src/rxode2-lincmt-carry-jump", compile = FALSE, quiet = TRUE)
  devtools::load_all("~/src/nlmixr2est-lincmt-speed", helpers = FALSE, quiet = TRUE)
})
rxode2::setRxThreads(1L)
nRep <- as.integer(Sys.getenv("REPS", "3"))

## ---- one-compartment oral model + data (same shape as the vs-ode bench) --
set.seed(1002003)
simTimes <- c(0.25, 0.75, 1.5, 3, 5.5, 8, 12.5, 17, 24, 32)
nSub <- 40L
eta <- matrix(rnorm(nSub * 3L, 0, 0.3), nSub, 3L,
              dimnames = list(NULL, c("ka", "cl", "v")))
simMod <- rxode2::rxode2(
  "cp = central/v; d/dt(depot) = -ka*depot; d/dt(central) = ka*depot - cl/v*central")
pars <- data.frame(ka = 1.2 * exp(eta[, "ka"]), cl = 4 * exp(eta[, "cl"]),
                   v = 30 * exp(eta[, "v"]))
ev <- rxode2::et(amt = 100, cmt = "depot") |> rxode2::et(simTimes)
sim <- rxode2::rxSolve(simMod, pars, ev, cores = 1L, addDosing = FALSE)
set.seed(2003004)
simDat <- data.frame(ID = rep(seq_len(nSub), each = length(simTimes)),
                     TIME = sim$time,
                     DV = sim$cp * (1 + rnorm(nrow(sim), 0, 0.15)),
                     AMT = 0, EVID = 0, CMT = "central")
doseRows <- data.frame(ID = seq_len(nSub), TIME = 0, DV = NA_real_,
                       AMT = 100, EVID = 1, CMT = "depot")
dat <- rbind(doseRows, simDat)
dat <- dat[order(dat$ID, dat$TIME, -dat$EVID), ]

# nolint start: object_usage_linter. (rxode2 ini/model DSL assignments)
uiLin <- function() {
  ini({
    lka <- log(1.44); lcl <- log(4.8); lv <- log(36)
    eta.ka ~ 0.1; eta.cl ~ 0.1; eta.v ~ 0.1
    prop.sd <- 0.2
  })
  model({
    ka <- exp(lka) * exp(eta.ka)
    cl <- exp(lcl) * exp(eta.cl)
    v <- exp(lv) * exp(eta.v)
    cp <- linCmt()
    cp ~ prop(prop.sd)
  })
}
# nolint end

ctlPost <- nlmixr2est::foceiControl(
  calcTables = FALSE, print = 0L, covMethod = "", maxOuterIterations = 0L,
  rxControl = rxode2::rxControl(cores = 1L, linCmtSensType = "AD"))

fitOnce <- function() {
  t0 <- proc.time()[["elapsed"]]
  suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(uiLin, dat, est = "focei", control = ctlPost)))
  c(sec = proc.time()[["elapsed"]] - t0)
}

## ---- cold vs warm wall times --------------------------------------------
secs <- numeric(0)
for (r in seq_len(nRep + 1L)) secs[r] <- fitOnce()  # rep 1 = cold (compiles)
cat(sprintf("cold (rep1, incl. compile): %.3f s\nwarm reps: %s (median %.3f s)\n",
            secs[1], paste(sprintf("%.3f", secs[-1]), collapse = ", "),
            stats::median(secs[-1])))

## ---- Rprof one warm repeat ----------------------------------------------
prof <- tempfile(fileext = ".out")
Rprof(prof, interval = 0.002, line.profiling = TRUE)
invisible(fitOnce())
Rprof(NULL)
s <- summaryRprof(prof)
tot <- s$sampling.time
cat(sprintf("\nprofiled warm fit: %.3f s sampled\n", tot))
cat("\n-- by.total, top 45 --\n")
print(utils::head(s$by.total, 45))
cat("\n-- by.self, top 30 --\n")
print(utils::head(s$by.self, 30))

## Stage aggregation: attribute by.total of DISJOINT anchor functions.
anchor <- c(
  uiBuild      = "rxode2::rxode2|rxUiCompress|rxUiDecompress",
  symengineGen = "\\brxS\\b|rxFromSE|rxOptExpr|rxNorm|\\.rxFinalizeInner|sympy|symengine",
  modelCompile = "rxModelVars|rxCompile|dynLoad|rxDll",
  etTransData  = "etTrans",
  foceiSetup   = "foceiSetup",
  innerEval    = "foceiInner|\\.Call.*Inner|foceiOuterF",
  tableTear    = "focei\\.theta|nlmixr2CreateOutputFromUi|\\.foceiFitInternal")
cat("\n-- anchor by.total seconds (overlapping; interpret with the tables above) --\n")
rn <- rownames(s$by.total)
for (a in names(anchor)) {
  hit <- grepl(anchor[[a]], rn)
  if (any(hit)) {
    cat(sprintf("%-14s max(by.total)=%6.3f s  fns: %s\n", a,
                max(s$by.total$total.time[hit]),
                paste(utils::head(gsub("\"", "", rn[hit]), 3), collapse = ", ")))
  } else cat(sprintf("%-14s (no samples)\n", a))
}
