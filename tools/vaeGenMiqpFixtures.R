# tools/vaeGenMiqpFixtures.R
# Generate the golden reference for the VAE covariate-selection kernel
# (vaeBestSubset_ in src/inner.cpp). The kernel finds the exact minimizer of the
# L0/BIC objective  score(S) = RSS_S / omega + penalty * |S|  that the reference
# Python MIQP (cvxpy + GUROBI, penalty = ln N) solves. For a moderate number of
# covariates that optimum is exactly reproducible by brute-force enumeration in
# R -- the SAME objective -- so the brute force is an authoritative, solver- and
# Python-independent oracle. The fixture pins several such cases (varying zDim,
# omega, penalty, fixed rows, empty/heavy supports) plus one large sparse case
# used to regress the nCov==32 overflow bug.
#
# Run:  Rscript tools/vaeGenMiqpFixtures.R
# Out:  tests/testthat/baselines/vae-miqp-golden.rds

set.seed(20260718L)

## Brute-force exact best subset for one response (independent of the C++ kernel).
bruteBestSubset <- function(y, X, omega, penalty) {
  nCov <- ncol(X)
  best <- list(score = Inf, sel = integer(0), coef = 0)
  for (m in 0:(2^nCov - 1)) {
    sel <- which(bitwAnd(m, bitwShiftL(1L, 0:(nCov - 1))) > 0)
    Xs <- cbind(1, X[, sel, drop = FALSE])
    coef <- qr.solve(Xs, y)
    rss <- sum((y - Xs %*% coef)^2)
    score <- rss / omega + penalty * length(sel)
    if (score < best$score - 1e-9) best <- list(score = score, sel = sel, coef = as.numeric(coef))
  }
  best
}

## Brute-force the whole [N, zDim] problem into intercept/beta/selected matrices.
bruteCase <- function(mu, X, omega, isFree, penalty) {
  zDim <- ncol(mu); nCov <- ncol(X)
  intercept <- numeric(zDim)
  beta <- matrix(0, zDim, nCov)
  selected <- matrix(0L, zDim, nCov)
  for (k in seq_len(zDim)) {
    if (isFree[k]) next
    b <- bruteBestSubset(mu[, k], X, omega[k], penalty)
    intercept[k] <- b$coef[1]
    if (length(b$sel) > 0) {
      beta[k, b$sel] <- b$coef[-1]
      selected[k, b$sel] <- 1L
    }
  }
  list(intercept = intercept, beta = beta, selected = selected)
}

makeCase <- function(name, N, nCov, zDim, omega, isFree, penalty, trueBeta, sd, seed) {
  set.seed(seed)
  X <- matrix(rnorm(N * nCov), N, nCov)
  mu <- matrix(0, N, zDim)
  for (k in seq_len(zDim)) {
    tb <- trueBeta[[k]]
    yk <- tb$intercept + (if (length(tb$sel)) X[, tb$sel, drop = FALSE] %*% tb$coef else 0)
    mu[, k] <- yk + rnorm(N, sd = sd)
  }
  list(name = name, N = N, nCov = nCov, zDim = zDim, penalty = penalty,
       inputs = list(mu = mu, covMat = X, omega = omega, isFree = isFree))
}

cases <- list()

## Case A: single eta, sparse true support {2, 5} out of 6.
cases$A <- makeCase("single-eta-6cov", N = 120L, nCov = 6L, zDim = 1L,
                    omega = 0.4, isFree = FALSE, penalty = log(120),
                    trueBeta = list(list(intercept = 1.0, sel = c(2L, 5L), coef = c(2.5, -1.8))),
                    sd = 0.5, seed = 101L)

## Case B: three etas, one free (skipped), different supports, 10 covariates.
cases$B <- makeCase("three-eta-10cov", N = 150L, nCov = 10L, zDim = 3L,
                    omega = c(0.3, 0.6, 0.2), isFree = c(FALSE, TRUE, FALSE),
                    penalty = log(150),
                    trueBeta = list(
                      list(intercept = 0.5, sel = c(1L, 7L), coef = c(2.0, 1.3)),
                      list(intercept = 0.0, sel = integer(0), coef = numeric(0)),
                      list(intercept = -0.7, sel = c(3L, 4L, 9L), coef = c(-1.5, 0.9, 2.2))),
                    sd = 0.4, seed = 202L)

## Case C: 14 covariates (2^14 subsets -- a real branch-and-bound stress test).
cases$C <- makeCase("single-eta-14cov", N = 200L, nCov = 14L, zDim = 1L,
                    omega = 0.35, isFree = FALSE, penalty = log(200),
                    trueBeta = list(list(intercept = 0.2, sel = c(4L, 8L, 12L),
                                         coef = c(2.4, -1.6, 1.1))),
                    sd = 0.45, seed = 303L)

## Case D: heavy penalty forces the empty model (no covariate should survive).
cases$D <- makeCase("empty-heavy-penalty", N = 100L, nCov = 8L, zDim = 1L,
                    omega = 0.5, isFree = FALSE, penalty = 1e6,
                    trueBeta = list(list(intercept = 1.0, sel = c(1L, 2L), coef = c(2.0, -2.0))),
                    sd = 0.5, seed = 404L)

## Compute the authoritative brute-force answer for the enumerable cases.
for (nm in names(cases)) {
  ci <- cases[[nm]]
  cases[[nm]]$expected <- bruteCase(ci$inputs$mu, ci$inputs$covMat, ci$inputs$omega,
                                    rep_len(ci$inputs$isFree, ci$zDim), ci$penalty)
  cat(sprintf("%-20s nCov=%2d zDim=%d selected: %s\n", ci$name, ci$nCov, ci$zDim,
              paste(apply(cases[[nm]]$expected$selected, 1,
                          function(r) paste(which(r == 1), collapse = ",")),
                    collapse = " | ")))
}

## Case E: nCov = 32 sparse strong signal -- brute force is infeasible, so the
## reference is the KNOWN generating support {3, 10, 25}. This regresses the
## `1u << 32` overflow (previously selected nothing) and proves scalability. The
## expected matrices are built directly from the ground truth; the test asserts
## the kernel recovers this support and runs quickly.
set.seed(505L)
N <- 200L; nCov <- 32L; selTrue <- c(3L, 10L, 25L); coefTrue <- c(2.0, -1.5, 1.2)
X32 <- matrix(rnorm(N * nCov), N, nCov)
y32 <- 0.5 + X32[, selTrue, drop = FALSE] %*% coefTrue + rnorm(N, sd = 0.4)
sel32 <- matrix(0L, 1L, nCov); sel32[1, selTrue] <- 1L
cases$E <- list(name = "sparse-32cov", N = N, nCov = nCov, zDim = 1L, penalty = log(N),
                inputs = list(mu = matrix(y32, ncol = 1), covMat = X32,
                              omega = 0.3, isFree = FALSE),
                expected = list(selected = sel32),  # support only (values not pinned)
                selectionOnly = TRUE)
cat(sprintf("%-20s nCov=%2d zDim=%d selected: %s (ground truth)\n",
            cases$E$name, cases$E$nCov, cases$E$zDim, paste(selTrue, collapse = ",")))

outDir <- file.path("tests", "testthat", "baselines")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
saveRDS(cases, file.path(outDir, "vae-miqp-golden.rds"), version = 2L)
cat("wrote", file.path(outDir, "vae-miqp-golden.rds"), "\n")
