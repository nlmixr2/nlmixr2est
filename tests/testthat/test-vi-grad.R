## FD gradient check for the ADVI outer likelihood + gradient (adviElboGrad_).
## At a FIXED eps draw the C++ analytic ELBO gradient (variational mu/omega,
## population thetas via the mu-ref identity + theta sensitivities, and the
## population between-subject log-variances) must match central finite
## differences of the returned ELBO.  Mean-field, point-estimate.

nmTest({
  test_that("adviElboGrad_ analytic gradient matches finite differences (mean-field)", {
    ## model with a mu-referenced structural theta (lka<->eta.ka), a non-mu
    ## structural theta (lcl, lv have no eta), and a residual sigma (add.err) --
    ## exercises all three population-gradient routes.
    mod <- function() {
      ini({ lka <- log(1.5); lcl <- log(0.09); lv <- log(32)
        eta.ka ~ 0.3; add.err <- 0.6 })
      model({ ka <- exp(lka + eta.ka); cl <- exp(lcl); v <- exp(lv)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - cl / v * central
        cp <- central / v; cp ~ add(add.err) })
    }
    ui <- rxode2::assertRxUi(mod)
    ## moderate ODE tolerance: tight enough for an accurate FD reference, but not
    ## so tight the 6-state sensitivity system trips the bad-solve tol relaxation.
    ctl <- emviControl(rxControl = rxode2::rxControl(atol = 1e-8, rtol = 1e-8))
    dat <- nlmixr2data::theo_sd
    N <- length(unique(dat$ID)); neta <- 1L

    ## classification -> mu-ref map (1-based ntheta per eta, NA if free)
    prep <- .adviDataPrep(ui, dat)
    muRefIdx <- prep$muRefThetaIdx
    muRefIdx[is.na(muRefIdx)] <- NA_integer_
    ntheta <- prep$ntheta

    ## fixed variational + population state and a fixed eps draw
    .testSeed(11)
    mu <- matrix(rnorm(N * neta, 0, 0.2), N, neta)
    omega <- matrix(rnorm(N * neta, -0.7, 0.1), N, neta)     # log-sd
    theta <- prep$theta                                      # natural-scale thetas
    logPopOmega <- log(prep$omega)                           # log population variances
    eps <- matrix(rnorm(N * neta), N, neta)

    .adviInnerSetup(ui, dat, mu, ctl)
    on.exit(.adviInnerFree(), add = TRUE)

    elboFun <- function(mu, omega, theta, logPopOmega) {
      adviElboGrad_(mu, omega, theta, logPopOmega, eps, as.integer(muRefIdx))$elbo
    }
    a <- adviElboGrad_(mu, omega, theta, logPopOmega, eps, as.integer(muRefIdx))
    ## central FD; the analytic gradient is exact, the FD is the noisy reference
    ## (ODE solve accuracy ~1e-7), so compare with a relative tolerance.
    h <- 1e-4
    fd1 <- function(f, x, i) {
      xp <- x; xm <- x; xp[i] <- x[i] + h; xm[i] <- x[i] - h
      (f(xp) - f(xm)) / (2 * h)
    }
    relOk <- function(a, fd, tol = 2e-3) abs(a - fd) < tol * (1 + abs(a))

    ## gMu / gOmega : a few subject/eta entries
    for (i in c(1L, 5L, 9L)) {
      gfd <- fd1(function(z) { m <- mu; m[i, 1] <- z[1]
        elboFun(m, omega, theta, logPopOmega) }, mu[i, 1, drop = FALSE], 1L)
      expect_true(relOk(a$gMu[i, 1], gfd))
      gfd <- fd1(function(z) { o <- omega; o[i, 1] <- z[1]
        elboFun(mu, o, theta, logPopOmega) }, omega[i, 1, drop = FALSE], 1L)
      expect_true(relOk(a$gOmega[i, 1], gfd))
    }

    ## gTheta : mu-ref (lka=1), non-mu struct (lcl=2, lv=3), sigma (add.err=4)
    for (p in seq_len(ntheta)) {
      gfd <- fd1(function(z) elboFun(mu, omega, z, logPopOmega), theta, p)
      expect_true(relOk(a$gTheta[p], gfd))
    }

    ## gPopLogOmega : population between-subject log-variance (no ODE dependence)
    for (k in seq_len(neta)) {
      gfd <- fd1(function(z) elboFun(mu, omega, theta, z), logPopOmega, k)
      expect_true(relOk(a$gPopLogOmega[k], gfd))
    }
  })

  test_that("adviElboGrad_ gTheta is right for an ESTIMATED boxCox lambda (#949)", {
    ## lambda is a residual-error theta that nevertheless reaches the MEAN: the
    ## ELBO depends on it through both sides of err = h(y; lambda) - h(f; lambda).
    ## Its theta-sensitivity column used to be identically zero, so the advi outer
    ## gradient had no lambda direction at all.  This is the advi-side counterpart
    ## of the C-API check in test-focei-ptrs.R -- adviElboGradCore reaches the same
    ## impThetaScore/impThetaAccumOne path the impmap M-step uses.
    mod <- function() {
      ini({ lka <- log(1.5); lcl <- log(0.09); lv <- log(32)
        eta.ka ~ 0.3; add.err <- 0.6; lambda <- c(-2, 0.8, 2) })
      model({ ka <- exp(lka + eta.ka); cl <- exp(lcl); v <- exp(lv)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - cl / v * central
        cp <- central / v; cp ~ add(add.err) + boxCox(lambda) })
    }
    ui <- rxode2::assertRxUi(mod)
    ctl <- emviControl(rxControl = rxode2::rxControl(atol = 1e-8, rtol = 1e-8))
    ## a boxCox endpoint floors a non-positive DV, which puts a kink in the
    ## objective that no finite difference can resolve -- theo_sd carries 9 such
    ## rows (DV == 0 pre-dose), really a censoring case.  Drop them for the check.
    dat <- nlmixr2data::theo_sd
    dat <- dat[!(dat$EVID == 0 & dat$DV <= 0), ]
    prep <- .adviDataPrep(ui, dat)
    N <- prep$N; neta <- 1L
    muRefIdx <- as.integer(prep$muRefThetaIdx)
    ntheta <- prep$ntheta
    ## which theta is lambda (do not hard-code the ini order)
    .th <- ui$iniDf[!is.na(ui$iniDf$ntheta), ]
    lamIdx <- .th$ntheta[.th$name == "lambda"]
    expect_length(lamIdx, 1L)

    .testSeed(11)
    mu <- matrix(rnorm(N * neta, 0, 0.2), N, neta)
    omega <- matrix(rnorm(N * neta, -0.7, 0.1), N, neta)
    theta <- prep$theta
    logPopOmega <- log(prep$omega)
    eps <- matrix(rnorm(N * neta), N, neta)

    .adviInnerSetup(ui, dat, mu, ctl)
    on.exit(.adviInnerFree(), add = TRUE)

    ## the DV-side term is actually WIRED, not merely giving the right number:
    ## lambdaOffset is >= 0 only when the model emits d(lambda)/d(theta)
    .info <- adviThetaSensInfo_()
    expect_true(.info$nSens > 0L)
    expect_true(.info$lambdaOffset >= 0L)

    elboFun <- function(theta)
      adviElboGrad_(mu, omega, theta, logPopOmega, eps, muRefIdx)$elbo
    a <- adviElboGrad_(mu, omega, theta, logPopOmega, eps, muRefIdx)
    h <- 1e-4
    relOk <- function(x, fd, tol = 2e-3) abs(x - fd) < tol * (1 + abs(x))
    for (p in seq_len(ntheta)) {
      xp <- theta; xp[p] <- theta[p] + h
      xm <- theta; xm[p] <- theta[p] - h
      gfd <- (elboFun(xp) - elboFun(xm)) / (2 * h)
      expect_true(relOk(a$gTheta[p], gfd),
                  info = paste0("theta ", p, ": analytic ", a$gTheta[p],
                                " fd ", gfd))
    }
    ## and the historical failure -- a lambda direction of exactly zero -- is gone
    expect_true(abs(a$gTheta[lamIdx]) > 1e-3)
  })

  test_that("adviElboGradFR_ full-rank gradient matches finite differences", {
    ## two correlated etas so the off-diagonal L entries are exercised
    mod <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; eta.cl ~ 0.3; add.err <- 0.6 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - cl / v * central
        cp <- central / v; cp ~ add(add.err) })
    }
    ui <- rxode2::assertRxUi(mod)
    ctl <- emviControl(rxControl = rxode2::rxControl(atol = 1e-8, rtol = 1e-8))
    dat <- nlmixr2data::theo_sd
    prep <- .adviDataPrep(ui, dat)
    N <- prep$N; neta <- 2L; nL <- neta * (neta + 1L) / 2L
    muRefIdx <- as.integer(prep$muRefThetaIdx); ntheta <- prep$ntheta

    .testSeed(2)
    mu <- matrix(rnorm(N * neta, 0, 0.2), N, neta)
    Lp <- matrix(rnorm(N * nL, 0, 0.05), N, nL)
    for (k in seq_len(neta)) Lp[, k * (k + 1L) / 2L] <- 0.4   # positive diagonal
    theta <- prep$theta; logPopOmega <- log(prep$omega); eps <- matrix(rnorm(N * neta), N, neta)

    .adviInnerSetup(ui, dat, mu, ctl)
    on.exit(.adviInnerFree(), add = TRUE)
    elbo <- function(mu, Lp, theta, lpo)
      adviElboGradFR_(mu, Lp, theta, lpo, eps, muRefIdx)$elbo
    a <- adviElboGradFR_(mu, Lp, theta, logPopOmega, eps, muRefIdx)
    h <- 1e-4; relOk <- function(x, fd) abs(x - fd) < 2e-3 * (1 + abs(x))

    ## gL : diagonal + off-diagonal entries for a couple of subjects
    for (s in c(1L, 5L, 9L)) for (j in seq_len(nL)) {
      gfd <- (elbo(mu, `[<-`(Lp, s, j, Lp[s, j] + h), theta, logPopOmega) -
              elbo(mu, `[<-`(Lp, s, j, Lp[s, j] - h), theta, logPopOmega)) / (2 * h)
      expect_true(relOk(a$gL[s, j], gfd))
    }
    ## gMu / gTheta / gPopLogOmega
    for (s in c(1L, 5L)) {
      gfd <- (elbo(`[<-`(mu, s, 1, mu[s, 1] + h), Lp, theta, logPopOmega) -
              elbo(`[<-`(mu, s, 1, mu[s, 1] - h), Lp, theta, logPopOmega)) / (2 * h)
      expect_true(relOk(a$gMu[s, 1], gfd))
    }
    for (p in seq_len(ntheta)) {
      gfd <- (elbo(mu, Lp, `[<-`(theta, p, theta[p] + h), logPopOmega) -
              elbo(mu, Lp, `[<-`(theta, p, theta[p] - h), logPopOmega)) / (2 * h)
      expect_true(relOk(a$gTheta[p], gfd))
    }
  })

  test_that("adviElboGrad_ gradient matches FD with a CORRELATED population omega", {
    ## With a declared eta block the population prior is the full quadratic form
    ## plus logdet(Omega), and gPopLogOmega becomes
    ##   0.5 w_k (Omega^-1 eta)_k^2 - 0.5 w_k (Omega^-1)_kk.
    ## FD-check that the returned ELBO and gradient stay mutually consistent --
    ## the diagonal tests above cannot see an error in either term.
    mod <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka + eta.cl ~ c(0.6,
                            0.1, 0.3)
        add.err <- 0.6 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) = -ka * depot; d/dt(central) = ka * depot - cl / v * central
        cp <- central / v; cp ~ add(add.err) })
    }
    ui <- rxode2::assertRxUi(mod)
    ctl <- emviControl(rxControl = rxode2::rxControl(atol = 1e-8, rtol = 1e-8))
    dat <- nlmixr2data::theo_sd
    prep <- .adviDataPrep(ui, dat)
    expect_true(.omegaHasOffDiag(prep$omegaMat))       # the block reached the prep
    N <- prep$N; neta <- 2L
    muRefIdx <- as.integer(prep$muRefThetaIdx); ntheta <- prep$ntheta

    .testSeed(7)
    mu <- matrix(rnorm(N * neta, 0, 0.2), N, neta)
    omega <- matrix(rnorm(N * neta, -0.7, 0.1), N, neta)     # log-sd
    theta <- prep$theta
    logPopOmega <- log(diag(prep$omegaMat))
    eps <- matrix(rnorm(N * neta), N, neta)

    .adviInnerSetup(ui, dat, mu, ctl)
    on.exit(.adviInnerFree(), add = TRUE)
    elboFun <- function(mu, omega, theta, lpo)
      adviElboGrad_(mu, omega, theta, lpo, eps, muRefIdx)$elbo
    a <- adviElboGrad_(mu, omega, theta, logPopOmega, eps, muRefIdx)
    h <- 1e-4
    fd1 <- function(f, x, i) {
      xp <- x; xm <- x; xp[i] <- x[i] + h; xm[i] <- x[i] - h
      (f(xp) - f(xm)) / (2 * h)
    }
    relOk <- function(x, fd, tol = 2e-3) abs(x - fd) < tol * (1 + abs(x))

    ## the population log-variance gradient under a correlated Omega
    for (k in seq_len(neta)) {
      gfd <- fd1(function(z) elboFun(mu, omega, theta, z), logPopOmega, k)
      expect_true(relOk(a$gPopLogOmega[k], gfd))
    }
    ## and the variational / theta gradients, which now see Omega^-1
    for (i in c(1L, 5L)) {
      gfd <- fd1(function(z) { m <- mu; m[i, 1] <- z[1]
        elboFun(m, omega, theta, logPopOmega) }, mu[i, 1, drop = FALSE], 1L)
      expect_true(relOk(a$gMu[i, 1], gfd))
    }
    for (p in seq_len(ntheta)) {
      gfd <- fd1(function(z) elboFun(mu, omega, z, logPopOmega), theta, p)
      expect_true(relOk(a$gTheta[p], gfd))
    }
  })
})
