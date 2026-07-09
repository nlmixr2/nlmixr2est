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
