# Measured trade-offs for the QRPEM options added in the Phoenix-parity work
# (developer script; not run on CI/CRAN).
#
#   NOT_CRAN=true Rscript design/qrpem/qrpem-options-bench.R
#
# Deliberately library(), NOT devtools::load_all(): load_all builds at -O0, so
# any timing it produced would not describe the shipping build, and the whole
# point of this script is to decide DEFAULTS for the shipping build.  Install
# first (R CMD INSTALL . or devtools::install()).
#
# Method, matching the convention the existing "Measured trade-off" blocks in
# R/impmap.R use: for each fixture, fit at a production isample over several
# seeds and score against a reference computed at a much larger isample.
# Reported per setting: theta RMSE, Omega RMSE, max Pareto k-hat, how many
# subjects exceed 0.7, mean effective-sample fraction, iterations, wall time.

library(nlmixr2)
library(nlmixr2est)

## ---- fixtures --------------------------------------------------------------
# 1 / 3 ETAs on theophylline, and an 8-ETA fixture: the eta dimension is what
# the scrambling hypothesis is actually about (a Cranley-Patterson shift is
# supposed to help least where the Sobol sequence degrades first), so a
# single-fixture answer would not settle it.
oneEta <- function() {
  ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; add.sd <- 0.7 })
  model({
    ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
    linCmt() ~ add(add.sd)
  })
}
threeEta <- function() {
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
# a general-likelihood (non-transformably-normal) endpoint: where the posteriors
# actually go non-Gaussian, and so where a laplace/mixture proposal has a job
llEta <- function() {
  ini({
    tka <- 0.45; tcl <- 1; tv <- 3.45
    eta.ka ~ 0.6; eta.cl ~ 0.3
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
    cp <- linCmt()
    ll(cp) ~ -0.5 * log(2 * pi) - log(add.sd) - 0.5 * ((DV - cp) / add.sd)^2
  })
}

fixtures <- list(eta1 = oneEta, eta3 = threeEta, ll2 = llEta)
dat <- nlmixr2data::theo_sd

## ---- settings under test ---------------------------------------------------
# Each is (label, extra impmapControl args).  "base" is the shipping default and
# every other row is scored against the same reference, so the comparison is
# tuned-against-tuned rather than one option against an unrelated baseline.
settings <- list(
  base         = list(),
  scrambleOwen = list(qr = TRUE, qrScramble = "owen"),
  scrambleLms  = list(qr = TRUE, qrScramble = "lms"),
  qrOnly       = list(qr = TRUE),
  laplace      = list(proposal = "laplace"),
  mixture      = list(proposal = "mixture"),
  burn5freeze  = list(nBurn = 5L, burnFreezeOmega = TRUE),
  mapIter3     = list(mapIter = 3L)
)

SEEDS   <- 1:8
NPROD   <- 300L      # production sample count
NREF    <- 8000L     # reference sample count
NITER   <- 60L

fitOne <- function(gen, seed, isample, nIter, extra = list()) {
  ctl <- do.call(impmapControl,
                 c(list(print = 0L, nIter = nIter, isample = as.integer(isample),
                        impSeed = as.integer(seed), covMethod = "",
                        calcTables = FALSE), extra))
  t0 <- proc.time()[["elapsed"]]
  f <- try(suppressWarnings(suppressMessages(nlmixr2(gen, dat, "impmap", ctl))),
           silent = TRUE)
  if (inherits(f, "try-error")) return(NULL)
  list(theta = fixef(f),
       omega = diag(f$omega),
       khat  = f$env$impPsisK,
       neff  = mean(f$env$impNeffFrac),
       iter  = f$env$impIter,
       secs  = proc.time()[["elapsed"]] - t0)
}

rmse <- function(x, ref) sqrt(mean((x - ref)^2))

score <- function(runs, refTheta, refOmega) {
  runs <- Filter(Negate(is.null), runs)
  if (!length(runs)) return(NULL)
  data.frame(
    n         = length(runs),
    thetaRMSE = mean(vapply(runs, function(r) rmse(r$theta, refTheta), numeric(1))),
    omegaRMSE = mean(vapply(runs, function(r) rmse(r$omega, refOmega), numeric(1))),
    maxKhat   = max(vapply(runs, function(r) max(r$khat, na.rm = TRUE), numeric(1))),
    nBadKhat  = mean(vapply(runs, function(r) sum(r$khat > 0.7, na.rm = TRUE), numeric(1))),
    neff      = mean(vapply(runs, function(r) r$neff, numeric(1))),
    iter      = mean(vapply(runs, function(r) r$iter, numeric(1))),
    secs      = mean(vapply(runs, function(r) r$secs, numeric(1)))
  )
}

## ---- run -------------------------------------------------------------------
out <- list()
for (fx in names(fixtures)) {
  gen <- fixtures[[fx]]
  cat("\n=== fixture:", fx, "===\n")
  cat("reference at isample =", NREF, "...\n")
  ref <- fitOne(gen, 1L, NREF, NITER)
  if (is.null(ref)) { cat("  reference FAILED, skipping fixture\n"); next }
  for (sname in names(settings)) {
    runs <- lapply(SEEDS, function(s)
      fitOne(gen, s, NPROD, NITER, settings[[sname]]))
    sc <- score(runs, ref$theta, ref$omega)
    if (is.null(sc)) { cat(sprintf("  %-13s ALL FAILED\n", sname)); next }
    sc$fixture <- fx; sc$setting <- sname
    out[[length(out) + 1L]] <- sc
    cat(sprintf("  %-13s thetaRMSE=%.5f omegaRMSE=%.5f maxK=%+.3f bad=%.2f neff=%.3f iter=%.1f %.1fs\n",
                sname, sc$thetaRMSE, sc$omegaRMSE, sc$maxKhat, sc$nBadKhat,
                sc$neff, sc$iter, sc$secs))
  }
}

res <- do.call(rbind, out)
res <- res[, c("fixture", "setting", "n", "thetaRMSE", "omegaRMSE", "maxKhat",
               "nBadKhat", "neff", "iter", "secs")]
print(res, row.names = FALSE)
saveRDS(res, "design/qrpem/qrpem-options-bench-results.rds")
cat("\nresults saved to design/qrpem/qrpem-options-bench-results.rds\n")
cat("\nA default moves only where a column wins here AND the loss elsewhere is\n")
cat("bounded -- weights with infinite variance are a correctness problem,\n")
cat("Monte-Carlo noise is a measurable cost.  See R/impmap.R's existing\n")
cat("'Measured trade-off' blocks for the format to paste these into.\n")
