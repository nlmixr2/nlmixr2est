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
  r <- vaeFoceLik(f0, J0, R, y, cens, limit, zLin, zPop, omega, off, 1, 0L, 30L, 1e-10, 1L)
  expect_lt(abs(r$objective - refObj(f0, J0, R, y, cens, limit, zLin, zPop, omega, off)), 1e-4)

  ## one M3 (BLQ) observation: obs 2 is below LOQ = its y value, cens = 1
  cens2 <- cens; cens2[2] <- 1L
  r2 <- vaeFoceLik(f0, J0, R, y, cens2, limit, zLin, zPop, omega, off, 1, 0L, 30L, 1e-10, 1L)
  expect_lt(abs(r2$objective - refObj(f0, J0, R, y, cens2, limit, zLin, zPop, omega, off)), 1e-4)
  ## the censored objective differs from treating the point as observed
  expect_gt(abs(r2$objective - r$objective), 1e-6)
})
