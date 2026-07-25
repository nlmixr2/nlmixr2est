## Wall-clock of one VAE covariate M-step selection, exact branch-and-bound vs
## L0Learn-proposed candidates + polish, at the neonatal case study's shape
## (N = 189 subjects).  The M-step runs one of these per latent dimension per
## training iteration, so the per-call cost below is multiplied by zDim * iters.
##
##   Rscript tools/benchVaeCovSelect.R
##
## Search size is reported in BITS of feasible-support space,
##   bits = sum over groups of log2(1 + columns in that group),
## which is the unit vaeControl(covSelectMaxExact=) is compared against.  One
## column per group is 1 bit (the historic candidate count); a covariate carrying
## both shape families is log2(3) = 1.585 bits.  Reporting bits is what makes the
## single-shape and multi-shape rows comparable, and is how the default
## covSelectMaxExact is calibrated.
##
## Not a test: it measures, it does not assert.

suppressMessages(devtools::load_all(".", helpers = FALSE, quiet = TRUE))
if (!requireNamespace("L0Learn", quietly = TRUE)) stop("benchmark needs L0Learn")

N <- 189L
set.seed(101)

## shapesPer = 1 reproduces the historic unconstrained sweep; shapesPer = 2 is the
## default multi-shape search (log + linear family per covariate).
run <- function(nGroups, shapesPer, bnbCapBits = 22) {
  nCov <- nGroups * shapesPer
  group <- rep(seq_len(nGroups), each = shapesPer)
  bits <- nGroups * log2(1 + shapesPer)
  X <- matrix(rnorm(N * nCov), N, nCov)
  ## columns within a group are alternate shapes of one covariate, so make them
  ## strongly correlated -- that is what the exclusion constraint has to arbitrate
  if (shapesPer > 1L) {
    for (g in seq_len(nGroups)) {
      .c <- which(group == g)
      for (.j in .c[-1]) X[, .j] <- X[, .c[1]] * 0.9 + rnorm(N, sd = 0.4)
    }
  }
  sel <- which(!duplicated(group) & group %in% sample.int(nGroups, min(3L, nGroups)))
  y <- as.numeric(0.5 + X[, sel, drop = FALSE] %*%
                    rep_len(c(1.5, -2, 1), length(sel)) + rnorm(N, sd = 1))
  omega <- 0.4
  penalty <- log(N)
  grpArg <- if (shapesPer > 1L) group else NULL

  t0 <- proc.time()[3]
  cand <- nlmixr2est:::.vaeL0Supports(X, y)
  got <- vaeScoreSupports_(y, X, omega, penalty, cand, polish = TRUE,
                           group = grpArg)
  tL <- proc.time()[3] - t0

  if (bits <= bnbCapBits) {
    t0 <- proc.time()[3]
    ref <- vaeBestSubset_(matrix(y, ncol = 1), X, omega, FALSE, penalty, "lifo",
                          grpArg)
    tB <- proc.time()[3] - t0
    same <- identical(which(got$selected == 1L), which(ref$selected[1, ] == 1L))
    cat(sprintf("%6d %7d %6d %8.2f %10.3f %10.4f %10.1f %8s\n",
                nGroups, shapesPer, nCov, bits, tB, tL, tB / tL, same))
  } else {
    cat(sprintf("%6d %7d %6d %8.2f %10s %10.4f %10s %8s\n",
                nGroups, shapesPer, nCov, bits, "skipped", tL, "-", "-"))
  }
}

cat(sprintf("%6s %7s %6s %8s %10s %10s %10s %8s\n",
            "groups", "shapes", "nCov", "bits", "bnb (s)", "l0 (s)",
            "speedup", "same"))
cat("-- one shape per covariate (historic search) --\n")
for (g in c(10L, 15L, 17L, 20L, 22L)) run(g, 1L)
cat("-- two shape families per covariate (new default) --\n")
for (g in c(6L, 8L, 10L, 11L, 12L, 13L)) run(g, 2L)
