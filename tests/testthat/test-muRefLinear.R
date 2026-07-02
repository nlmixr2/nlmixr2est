test_that(".muRefLin recovers known population theta/covariate coefficient", {
  # Pure R math prototype (Phase 0) -- no compiled code, no fitting.
  # Simulate subjects with a known population intercept/slope and
  # genuine random noise; with enough subjects OLS should recover both
  # within a reasonable tolerance.
  set.seed(42)
  n <- 300
  trueTheta <- 1.5
  trueBeta <- 0.5
  cov <- data.frame(logWT = rnorm(n, 0, 0.3))
  eta <- rnorm(n, 0, 0.2)
  phi <- trueTheta + trueBeta * cov$logWT + eta

  res <- nlmixr2est:::.muRefLin(phi, cov)
  expect_equal(res$theta, trueTheta, tolerance = 0.05)
  expect_equal(unname(res$coef["logWT"]), trueBeta, tolerance = 0.1)

  # round-trip identity: theta + coef*cov + eta must reproduce phi exactly
  # (this holds regardless of how well theta/coef happen to match the
  # simulation truth -- it's an algebraic property of OLS, not a
  # statistical one)
  reconstructed <- res$theta + res$coef["logWT"] * cov$logWT + res$eta
  expect_equal(unname(reconstructed), phi, tolerance = 1e-8)
})

test_that(".muRefLin respects a user-fixed covariate coefficient", {
  set.seed(43)
  n <- 100
  trueTheta <- 0.8
  fixedBeta <- 0.75
  cov <- data.frame(logWT = rnorm(n, 0, 0.3))
  eta <- rnorm(n, 0, 0.1)
  phi <- trueTheta + fixedBeta * cov$logWT + eta

  res <- nlmixr2est:::.muRefLin(phi, cov, fixedCoef = c(logWT = fixedBeta))
  # the fixed coefficient is returned unchanged, not re-estimated
  expect_equal(unname(res$coef["logWT"]), fixedBeta)
  expect_equal(res$theta, trueTheta, tolerance = 0.05)
  reconstructed <- res$theta + res$coef["logWT"] * cov$logWT + res$eta
  expect_equal(unname(reconstructed), phi, tolerance = 1e-8)
})

test_that(".muRefLin handles multiple covariates, some fixed and some free", {
  set.seed(44)
  n <- 200
  trueTheta <- 1.0
  trueBetaFree <- 0.4
  fixedBeta <- -0.3
  cov <- data.frame(logWT = rnorm(n, 0, 0.3), sexf = rbinom(n, 1, 0.5))
  eta <- rnorm(n, 0, 0.15)
  phi <- trueTheta + trueBetaFree * cov$logWT + fixedBeta * cov$sexf + eta

  res <- nlmixr2est:::.muRefLin(phi, cov, fixedCoef = c(sexf = fixedBeta))
  expect_equal(unname(res$coef["sexf"]), fixedBeta)
  expect_equal(unname(res$coef["logWT"]), trueBetaFree, tolerance = 0.1)
  expect_equal(res$theta, trueTheta, tolerance = 0.1)
})

test_that(".muRefIrls with equal weights matches .muRefLin exactly", {
  set.seed(45)
  n <- 80
  cov <- data.frame(logWT = rnorm(n, 0, 0.3))
  phi <- 1.2 + 0.6 * cov$logWT + rnorm(n, 0, 0.2)

  resLin <- nlmixr2est:::.muRefLin(phi, cov)
  resIrls <- nlmixr2est:::.muRefIrls(phi, cov, weights = rep(1, n))
  expect_equal(resIrls$theta, resLin$theta, tolerance = 1e-10)
  expect_equal(unname(resIrls$coef["logWT"]), unname(resLin$coef["logWT"]), tolerance = 1e-10)
  expect_equal(resIrls$eta, resLin$eta, tolerance = 1e-10)
})

test_that(".muRefIrls downweights high-variance subjects toward the low-variance estimate", {
  # Two-group heteroscedastic design: half the subjects have small residual
  # variance (precise phi), half have large (noisy phi). A properly
  # curvature-weighted regression should track the precise group more
  # closely than an unweighted OLS fit would.
  set.seed(46)
  n <- 400
  # covariate spans the same full range in both the precise and noisy
  # halves (interleaved, not split by covariate value) so the slope stays
  # identifiable even after the noisy half is downweighted -- a design
  # split by covariate *level* instead would make the slope itself
  # unidentifiable once one whole level is downweighted, which would test
  # the design, not the weighting.
  cov <- data.frame(logWT = rep(seq(-0.3, 0.3, length.out = n / 2), 2))
  trueTheta <- 1.0
  trueBeta <- 0.5
  precise <- seq_len(n) <= n / 2
  sdVec <- ifelse(precise, 0.05, 2)
  eta <- rnorm(n, 0, sdVec)
  phi <- trueTheta + trueBeta * cov$logWT + eta
  weights <- 1 / sdVec^2

  resLin <- nlmixr2est:::.muRefLin(phi, cov)
  resIrls <- nlmixr2est:::.muRefIrls(phi, cov, weights = weights)

  expect_equal(resIrls$theta, trueTheta, tolerance = 0.05)
  expect_equal(unname(resIrls$coef["logWT"]), trueBeta, tolerance = 0.1)
  # the correctly-weighted fit should be at least as close to the truth as
  # plain OLS on this sharply heteroscedastic design
  expect_true(abs(resIrls$theta - trueTheta) <= abs(resLin$theta - trueTheta) + 1e-8)
})
