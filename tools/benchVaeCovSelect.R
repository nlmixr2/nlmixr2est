## Wall-clock of one VAE covariate M-step selection, exact branch-and-bound vs
## L0Learn-proposed candidates + polish, at the neonatal case study's shape
## (N = 189 subjects).  The M-step runs one of these per latent dimension per
## training iteration, so the per-call cost below is multiplied by zDim * iters.
##
##   Rscript tools/benchVaeCovSelect.R
##
## Not a test: it measures, it does not assert.

suppressMessages(devtools::load_all(".", helpers = FALSE, quiet = TRUE))
if (!requireNamespace("L0Learn", quietly = TRUE)) stop("benchmark needs L0Learn")

N <- 189L
nCovs <- c(10L, 15L, 20L, 25L, 30L, 35L)
bnbCap <- 30L                     # above this the exact search is left unmeasured
set.seed(101)

cat(sprintf("%5s %10s %10s %10s %8s\n", "nCov", "bnb (s)", "l0 (s)", "speedup", "same"))
for (nCov in nCovs) {
  X <- matrix(rnorm(N * nCov), N, nCov)
  sel <- sort(sample.int(nCov, 3L))
  y <- as.numeric(0.5 + X[, sel, drop = FALSE] %*% c(1.5, -2, 1) + rnorm(N, sd = 1))
  omega <- 0.4
  penalty <- log(N)

  t0 <- proc.time()[3]
  cand <- nlmixr2est:::.vaeL0Supports(X, y)
  got <- vaeScoreSupports_(y, X, omega, penalty, cand, polish = TRUE)
  tL <- proc.time()[3] - t0

  if (nCov <= bnbCap) {
    t0 <- proc.time()[3]
    ref <- vaeBestSubset_(matrix(y, ncol = 1), X, omega, FALSE, penalty)
    tB <- proc.time()[3] - t0
    same <- identical(which(got$selected == 1L), which(ref$selected[1, ] == 1L))
    cat(sprintf("%5d %10.3f %10.4f %10.1f %8s\n", nCov, tB, tL, tB / tL, same))
  } else {
    cat(sprintf("%5d %10s %10.4f %10s %8s\n", nCov, "skipped", tL, "-", "-"))
  }
}
