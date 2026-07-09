# RPEM likelihood model (design/rpem/13, criteria C1.x in design/rpem/14).
# The model must return log p(Y_i | theta_i) for a supplied theta = mu + eta
# with no eta optimization.

.rpemOneCmt <- function() {
  ini({
    tka <- 0.45; tcl <- 1.0; tv <- 3.45; add.sd <- 0.7
    eta.ka ~ 0.6
  })
  model({
    ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

.rpemLlikData <- function() {
  ev <- rxode2::et(amt=100, cmt="depot")
  ev <- rxode2::et(ev, seq(0.5, 24, by=2.5))
  dat <- as.data.frame(ev)
  dat$DV <- 0
  .obs <- dat$evid == 0
  dat$DV[.obs] <- c(4, 7, 6, 5, 4, 3, 2.5, 2, 1.5, 1.2)[seq_len(sum(.obs))]
  dat
}

.rpemSolve <- function(m, th, etaVal, dat) {
  .pars <- stats::setNames(c(th, etaVal),
                           c("THETA[1]", "THETA[2]", "THETA[3]", "THETA[4]", "ETA[1]"))
  rxode2::rxSolve(m, params=.pars, events=dat, returnType="data.frame",
                  addDosing=FALSE, atol=1e-10, rtol=1e-10)
}

test_that("rpem likelihood model returns log p(Y|theta) for supplied eta (C1.x)", {
  skip_on_cran()
  ui <- rxode2::rxode2(.rpemOneCmt)
  m  <- ui$rpemRxModel$predOnly
  th <- c(0.45, 1.0, 3.45, 0.7)
  addSd <- th[4]
  dat <- .rpemLlikData()
  dvObs <- dat$DV[dat$evid == 0]

  # C1.1: -rx_pred_ equals the additive-error normal log-likelihood at cp
  s <- .rpemSolve(m, th, 0.1, dat)
  cp <- s$rx_pred_f_
  refll <- dnorm(dvObs, mean=cp, sd=addSd, log=TRUE)
  expect_equal(-s$rx_pred_, refll, tolerance=1e-8)

  # C1.2: eta=0 structural prediction reproduces the population solve
  s0 <- .rpemSolve(m, th, 0, dat)
  ref <- rxode2::rxode2({
    ka <- exp(0.45); cl <- exp(1.0); v <- exp(3.45)
    cp <- linCmt()
  })
  sref <- rxode2::rxSolve(ref, events=dat, returnType="data.frame",
                          addDosing=FALSE, atol=1e-10, rtol=1e-10)
  expect_equal(s0$rx_pred_f_, sref$cp, tolerance=1e-9)

  # C1.3: residual-only re-scoring (reuse cached cp) matches a full re-solve
  newSd <- 1.3
  llFromCachedCp <- dnorm(dvObs, mean=cp, sd=newSd, log=TRUE)
  sNew <- .rpemSolve(m, c(th[1:3], newSd), 0.1, dat)
  expect_equal(-sNew$rx_pred_, llFromCachedCp, tolerance=1e-8)
})
