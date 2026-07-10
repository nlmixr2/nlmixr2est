# RPEM C++ E-step solve (src/rpem.cpp): rpemSetup + rpemSolvePop must return
# per-subject -sum(rx_pred_) == -log p(Y_i | theta_i), matching the R rxSolve
# path (design/rpem/04, criteria C1.x / C8.1).

test_that("rpem C++ solve matches R likelihood per subject", {
  skip_on_cran()
  one.cmt <- function() {
    ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; add.sd <- 0.7; eta.ka ~ 0.6 })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
      cp <- linCmt(); cp ~ add(add.sd)
    })
  }
  ui <- rxode2::rxode2(one.cmt)
  m  <- ui$rpemRxModel$predOnly

  ev <- rxode2::et(amt=100, cmt="depot")
  ev <- rxode2::et(ev, seq(0.5, 24, by=2.5))
  dat <- as.data.frame(ev); dat$DV <- 0
  .obs <- dat$evid == 0
  dat$DV[.obs] <- c(4, 7, 6, 5, 4, 3, 2.5, 2, 1.5, 1.2)[seq_len(sum(.obs))]

  th <- c(0.45, 1.0, 3.45, 0.7); eta <- 0.1

  pars <- stats::setNames(c(th, eta),
                          c("THETA[1]", "THETA[2]", "THETA[3]", "THETA[4]", "ETA[1]"))
  sR <- rxode2::rxSolve(m, params=pars, events=dat, returnType="data.frame",
                        addDosing=FALSE, atol=1e-10, rtol=1e-10)
  sumPredR <- sum(sR$rx_pred_)

  e <- new.env()
  e$predOnly <- m
  e$rxControl <- rxode2::rxControl(atol=1e-10, rtol=1e-10)
  e$param <- stats::setNames(c(th, eta),
                             c("THETA[1]", "THETA[2]", "THETA[3]", "THETA[4]", "ETA[1]"))
  e$data  <- dat
  rpemSetup(e)
  on.exit(rpemFree(), add=TRUE)
  sumPredC <- rpemSolvePop(matrix(c(th, eta), nrow=1))

  expect_equal(sumPredC[1], sumPredR, tolerance=1e-6)
})

test_that("rpem C++ E-step K=1 log-sum-exp accumulation matches R (Eq 24/26)", {
  skip_on_cran()
  one.cmt <- function() {
    ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; add.sd <- 0.7; eta.ka ~ 0.6 })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
      cp <- linCmt(); cp ~ add(add.sd)
    })
  }
  ui <- rxode2::rxode2(one.cmt)
  m  <- ui$rpemRxModel$predOnly

  ev <- rxode2::et(amt=100, cmt="depot")
  ev <- rxode2::et(ev, seq(0.5, 24, by=2.5))
  dat <- as.data.frame(ev); dat$DV <- 0
  .obs <- dat$evid == 0
  dat$DV[.obs] <- c(4, 7, 6, 5, 4, 3, 2.5, 2, 1.5, 1.2)[seq_len(sum(.obs))]

  th <- c(0.45, 1.0, 3.45, 0.7); om <- 0.6
  set.seed(42); nG <- 200
  etas <- rnorm(nG, 0, sqrt(om))
  parBig <- cbind(matrix(rep(th, each=nG), nrow=nG), etas) # nG x 5 (THETA1..4, ETA1)

  e <- new.env()
  e$predOnly <- m
  e$rxControl <- rxode2::rxControl(atol=1e-10, rtol=1e-10)
  e$param <- stats::setNames(c(th, 0),
                             c("THETA[1]", "THETA[2]", "THETA[3]", "THETA[4]", "ETA[1]"))
  e$data  <- dat
  rpemSetup(e)
  on.exit(rpemFree(), add=TRUE)
  res <- rpemEstepK1(parBig, nG)

  # R reference: same draws, log p per draw, stable log-sum-exp - log(nG)
  logp <- vapply(seq_len(nG), function(j) {
    pars <- stats::setNames(c(th, etas[j]),
                            c("THETA[1]", "THETA[2]", "THETA[3]", "THETA[4]", "ETA[1]"))
    s <- rxode2::rxSolve(m, params=pars, events=dat, returnType="data.frame",
                         addDosing=FALSE, atol=1e-10, rtol=1e-10)
    -sum(s$rx_pred_)
  }, numeric(1))
  mx <- max(logp)
  logniR <- mx + log(sum(exp(logp - mx))) - log(nG)

  expect_equal(res$lnL, logniR, tolerance=1e-8)
  expect_equal(res$logn[1], logniR, tolerance=1e-8)
})

test_that("rpem C++ E-step draws etas via threefry rxRmvn (D18) and is reproducible", {
  skip_on_cran()
  one.cmt <- function() {
    ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; add.sd <- 0.7; eta.ka ~ 0.6 })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
      cp <- linCmt(); cp ~ add(add.sd)
    })
  }
  ui <- rxode2::rxode2(one.cmt)
  m  <- ui$rpemRxModel$predOnly

  ev <- rxode2::et(amt=100, cmt="depot")
  ev <- rxode2::et(ev, seq(0.5, 24, by=2.5))
  dat <- as.data.frame(ev); dat$DV <- 0
  .obs <- dat$evid == 0
  dat$DV[.obs] <- c(4, 7, 6, 5, 4, 3, 2.5, 2, 1.5, 1.2)[seq_len(sum(.obs))]

  th <- c(0.45, 1.0, 3.45, 0.7); om <- 0.6; nG <- 200

  e <- new.env()
  e$predOnly <- m
  e$rxControl <- rxode2::rxControl(atol=1e-10, rtol=1e-10)
  e$param <- stats::setNames(c(th, 0),
                             c("THETA[1]", "THETA[2]", "THETA[3]", "THETA[4]", "ETA[1]"))
  e$data  <- dat
  rpemSetup(e)
  on.exit(rpemFree(), add=TRUE)

  base <- c(th, 0)         # ETA[1] slot placeholder, overwritten per draw
  etaIdx <- 4L             # 0-based position of ETA[1] in the param row
  omega <- matrix(om, 1, 1)

  rxode2::rxSetSeed(1234)
  etaMat <- rxode2::rxRmvn(nG, mu = 0, sigma = omega)   # single subject: nAll = nG
  res <- rpemEstepK1Draw(e, base, etaIdx, etaMat, nG, 1L, matrix(0, nrow(etaMat)/nG, ncol(etaMat)), rep(1, ncol(etaMat)), 1.0)

  # correctness: recompute lnL from the etas the C++ engine drew
  etas <- as.numeric(res$eta)
  logp <- vapply(seq_len(nG), function(j) {
    pars <- stats::setNames(c(th, etas[j]),
                            c("THETA[1]", "THETA[2]", "THETA[3]", "THETA[4]", "ETA[1]"))
    s <- rxode2::rxSolve(m, params=pars, events=dat, returnType="data.frame",
                         addDosing=FALSE, atol=1e-10, rtol=1e-10)
    -sum(s$rx_pred_)
  }, numeric(1))
  mx <- max(logp)
  lnLref <- mx + log(sum(exp(logp - mx))) - log(nG)
  expect_equal(res$lnL, lnLref, tolerance=1e-8)

  # reproducibility (C3.1): same threefry seed -> identical draws and lnL
  rxode2::rxSetSeed(1234)
  etaMat2 <- rxode2::rxRmvn(nG, mu = 0, sigma = omega)
  res2 <- rpemEstepK1Draw(e, base, etaIdx, etaMat2, nG, 1L, matrix(0, nrow(etaMat2)/nG, ncol(etaMat2)), rep(1, ncol(etaMat2)), 1.0)
  expect_equal(res2$eta, res$eta)
  expect_equal(res2$lnL, res$lnL)
})

test_that("rpem C++ E-step is correct and deterministic (multi-subject)", {
  skip_on_cran()
  one.cmt <- function() {
    ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; add.sd <- 0.7; eta.ka ~ 0.6 })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
      cp <- linCmt(); cp ~ add(add.sd)
    })
  }
  ui <- rxode2::rxode2(one.cmt)
  m  <- ui$rpemRxModel$predOnly

  nid <- 3L
  dat <- do.call(rbind, lapply(seq_len(nid), function(i) {
    ev <- rxode2::et(amt=100, cmt="depot")
    ev <- rxode2::et(ev, seq(0.5, 24, by=2.5))
    d <- as.data.frame(ev); d$id <- i; d$DV <- 0
    o <- d$evid == 0
    d$DV[o] <- (c(4, 7, 6, 5, 4, 3, 2.5, 2, 1.5, 1.2) + i)[seq_len(sum(o))]
    d
  }))

  th <- c(0.45, 1.0, 3.45, 0.7); om <- 0.6; nG <- 100
  nm <- c("THETA[1]", "THETA[2]", "THETA[3]", "THETA[4]", "ETA[1]")

  e <- new.env()
  e$predOnly <- m
  e$rxControl <- rxode2::rxControl(atol=1e-10, rtol=1e-10)
  e$param <- stats::setNames(c(th, 0), nm)
  e$data  <- dat
  rpemSetup(e)
  on.exit(rpemFree(), add=TRUE)

  base <- c(th, 0); etaIdx <- 4L; omega <- matrix(om, 1, 1)

  rxode2::rxSetSeed(7)
  etaMat <- rxode2::rxRmvn(nid * nG, mu = 0, sigma = omega)
  r1 <- rpemEstepK1Draw(e, base, etaIdx, etaMat, nG, 2L, matrix(0, nrow(etaMat)/nG, ncol(etaMat)), rep(1, ncol(etaMat)), 1.0)   # solve on 2 cores

  # correctness under threading: recompute each subject's logn from its etas
  etas <- as.numeric(r1$eta)
  lognRef <- vapply(seq_len(nid), function(i) {
    dsub <- dat[dat$id == i, ]
    idx <- ((i - 1L) * nG + 1L):(i * nG)
    logp <- vapply(idx, function(r) {
      pars <- stats::setNames(c(th, etas[r]), nm)
      s <- rxode2::rxSolve(m, params=pars, events=dsub, returnType="data.frame",
                           addDosing=FALSE, atol=1e-10, rtol=1e-10)
      -sum(s$rx_pred_)
    }, numeric(1))
    mx <- max(logp); mx + log(sum(exp(logp - mx))) - log(nG)
  }, numeric(1))
  expect_equal(r1$logn, lognRef, tolerance=1e-8)

  # determinism (C3.1): same seed + same cores -> identical
  rxode2::rxSetSeed(7)
  etaMat2 <- rxode2::rxRmvn(nid * nG, mu = 0, sigma = omega)
  r2 <- rpemEstepK1Draw(e, base, etaIdx, etaMat2, nG, 2L, matrix(0, nrow(etaMat2)/nG, ncol(etaMat2)), rep(1, ncol(etaMat2)), 1.0)
  expect_equal(r2$logn, r1$logn)
  expect_equal(r2$lnL, r1$lnL)
})
