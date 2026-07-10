## Phase 6 polish: the C++ FOCE-linearized marginal -2LL kernel (vaeFoceLik) with
## MAP-EBE re-optimization and M2/M3/M4 censoring. Validated against an
## independent R implementation of the same Laplace objective (with numeric
## per-obs derivatives), for an uncensored problem and one with an M3 (BLQ)
## observation -- confirming the censored contribution matches log Phi((LOQ-f)/sd).

test_that("vaeFoceLik matches an R Laplace reference (observed and M3-censored)", {
  ln2pi <- log(2 * pi)
  ## per-observation negative log-likelihood: observed Gaussian or M3 censored
  negll <- function(f, R, y, cens, lim) {
    sd <- sqrt(R)
    if (cens == 0 && !is.finite(lim)) return(0.5 * ((y - f)^2 / R + log(R) + ln2pi))
    -stats::pnorm(cens * (y - f) / sd, log.p = TRUE)   # M3
  }
  ## R reference: per-subject Newton to the FOCE mode of the linearized model,
  ## then the Laplace -2LL. Derivatives of negll by central differences.
  refObj <- function(f0, J0, R, y, cens, limit, zLin, zPop, omega, off) {
    N <- nrow(zLin); neta <- ncol(zLin); Oinv <- diag(1 / omega, neta)
    total <- 0
    for (i in seq_len(N)) {
      idx <- (off[i] + 1L):off[i + 1L]
      z <- zLin[i, ]
      dObs <- function(f, j) {
        h <- 1e-5
        g <- (negll(f + h, R[j], y[j], cens[j], limit[j]) - negll(f - h, R[j], y[j], cens[j], limit[j])) / (2 * h)
        hh <- (negll(f + h, R[j], y[j], cens[j], limit[j]) - 2 * negll(f, R[j], y[j], cens[j], limit[j]) +
                 negll(f - h, R[j], y[j], cens[j], limit[j])) / h^2
        c(g, hh)
      }
      for (it in 1:30) {
        grad <- Oinv %*% (z - zPop[i, ]); H <- Oinv
        for (j in idx) {
          Jj <- J0[j, ]; f <- f0[j] + sum(Jj * (z - zLin[i, ])); d <- dObs(f, j)
          grad <- grad + d[1] * Jj; H <- H + d[2] * outer(Jj, Jj)
        }
        step <- solve(H, grad); z <- z - as.numeric(step)
        if (sqrt(sum(step^2)) < 1e-10) break
      }
      ll <- 0; H <- Oinv
      for (j in idx) {
        Jj <- J0[j, ]; f <- f0[j] + sum(Jj * (z - zLin[i, ]))
        ll <- ll - negll(f, R[j], y[j], cens[j], limit[j])
        H <- H + dObs(f, j)[2] * outer(Jj, Jj)
      }
      ez <- z - zPop[i, ]
      total <- total - 2 * ll + sum(ez^2 / omega) + sum(log(omega)) + as.numeric(determinant(H, TRUE)$modulus)
    }
    total
  }

  set.seed(3)
  N <- 3L; neta <- 2L; ni <- 4L
  off <- as.integer(c(0, cumsum(rep(ni, N))))
  tot <- N * ni
  f0 <- runif(tot, 1, 5); J0 <- matrix(rnorm(tot * neta, 0, 0.5), tot, neta)
  R <- rep(0.2, tot); y <- f0 + rnorm(tot, 0, 0.3)
  zLin <- matrix(rnorm(N * neta, 0, 0.2), N, neta); zPop <- matrix(0, N, neta)
  omega <- c(0.3, 0.15)

  ## uncensored
  cens <- integer(tot); limit <- rep(NA_real_, tot)
  r <- vaeFoceLik(f0, J0, R, y, cens, limit, zLin, zPop, omega, off, 1, 0L, 30L, 1e-10, 1L, 1L, 1)
  expect_lt(abs(r$objective - refObj(f0, J0, R, y, cens, limit, zLin, zPop, omega, off)), 1e-4)

  ## one M3 (BLQ) observation: obs 2 is below LOQ = its y value, cens = 1
  cens2 <- cens; cens2[2] <- 1L
  r2 <- vaeFoceLik(f0, J0, R, y, cens2, limit, zLin, zPop, omega, off, 1, 0L, 30L, 1e-10, 1L, 1L, 1)
  expect_lt(abs(r2$objective - refObj(f0, J0, R, y, cens2, limit, zLin, zPop, omega, off)), 1e-4)
  ## the censored objective differs from treating the point as observed
  expect_gt(abs(r2$objective - r$objective), 1e-6)
})

test_that("vaeFoceLik combines mixture components as -2 logsumexp (inner.cpp)", {
  ln2pi <- log(2 * pi)
  negll <- function(f, R, y) 0.5 * ((y - f)^2 / R + log(R) + ln2pi)
  ## per-subject Laplace -2LL for one component (uncensored, quadratic -> exact)
  subjObj <- function(f0, J0, R, y, zLin, zPop, omega, idx) {
    neta <- length(zLin); Oinv <- diag(1 / omega, neta); z <- zLin
    for (it in 1:50) {
      grad <- Oinv %*% (z - zPop); H <- Oinv
      for (j in idx) {
        Jj <- J0[j, ]; f <- f0[j] + sum(Jj * (z - zLin))
        grad <- grad + (-(y[j] - f) / R[j]) * Jj; H <- H + (1 / R[j]) * outer(Jj, Jj)
      }
      step <- solve(H, grad); z <- z - as.numeric(step)
      if (sqrt(sum(step^2)) < 1e-12) break
    }
    ll <- 0; H <- Oinv
    for (j in idx) {
      Jj <- J0[j, ]; f <- f0[j] + sum(Jj * (z - zLin))
      ll <- ll - negll(f, R[j], y[j]); H <- H + (1 / R[j]) * outer(Jj, Jj)
    }
    ez <- z - zPop
    -2 * ll + sum(ez^2 / omega) + sum(log(omega)) + as.numeric(determinant(H, TRUE)$modulus)
  }

  set.seed(11)
  N <- 3L; neta <- 2L; ni <- 4L; nMix <- 2L
  off <- as.integer(c(0, cumsum(rep(ni, N)))); tot <- N * ni
  y <- runif(tot, 1, 5)
  ## component-major stacked data (component 2 predictions shifted)
  f0c <- list(runif(tot, 1, 5), runif(tot, 1, 5) + 0.4)
  J0c <- list(matrix(rnorm(tot * neta, 0, 0.5), tot, neta),
              matrix(rnorm(tot * neta, 0, 0.5), tot, neta))
  Rc <- list(rep(0.2, tot), rep(0.25, tot))
  zLc <- list(matrix(rnorm(N * neta, 0, 0.2), N, neta), matrix(rnorm(N * neta, 0, 0.2), N, neta))
  zPc <- list(matrix(0, N, neta), matrix(0.05, N, neta))
  omega <- c(0.3, 0.15); mixProb <- c(0.65, 0.35)

  f0 <- do.call(c, f0c); J0 <- do.call(rbind, J0c); R <- do.call(c, Rc)
  zLin <- do.call(rbind, zLc); zPop <- do.call(rbind, zPc)
  cens <- integer(tot); limit <- rep(NA_real_, tot)

  r <- vaeFoceLik(f0, J0, R, y, cens, limit, zLin, zPop, omega, off, 1, 0L, 30L, 1e-12,
                  1L, nMix, mixProb)

  ## R reference: per subject -2 logsumexp_m( log p_m - 0.5 obj_i^(m) )
  ref <- 0
  for (i in seq_len(N)) {
    idx <- (off[i] + 1L):off[i + 1L]
    llc <- vapply(seq_len(nMix), function(m) {
      log(mixProb[m]) - 0.5 * subjObj(f0c[[m]], J0c[[m]], Rc[[m]], y, zLc[[m]][i, ], zPc[[m]][i, ], omega, idx)
    }, numeric(1))
    mmax <- max(llc)
    ref <- ref - 2 * (mmax + log(sum(exp(llc - mmax))))
  }
  expect_lt(abs(r$objective - ref), 1e-4)

  ## nMix=1 with mixProb=1 equals the single-component objective (component 1)
  r1 <- vaeFoceLik(f0c[[1]], J0c[[1]], Rc[[1]], y, cens, limit, zLc[[1]], zPc[[1]],
                   omega, off, 1, 0L, 30L, 1e-12, 1L, 1L, 1)
  s1 <- sum(vapply(seq_len(N), function(i) {
    idx <- (off[i] + 1L):off[i + 1L]
    subjObj(f0c[[1]], J0c[[1]], Rc[[1]], y, zLc[[1]][i, ], zPc[[1]][i, ], omega, idx)
  }, numeric(1)))
  expect_lt(abs(r1$objective - s1), 1e-4)
  ## mixest reports a valid component per subject
  expect_true(all(r$mixest %in% seq_len(nMix)))
})
