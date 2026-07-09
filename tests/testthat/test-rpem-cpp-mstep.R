# RPEM C++ M-step (src/rpem.cpp, design/rpem/05).  The joint (subject,sample)
# Metropolis-Hastings walk over the stored E-step samples must reproduce the
# conjugate mu/Omega update (Eq 15-16), which for K=1 equals the self-normalized
# importance-sampling estimate over the same samples.

test_that("rpem C++ M-step MH matches importance-sampling conjugate update (K=1)", {
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

  th <- c(0.45, 1.0, 3.45, 0.7); om <- 0.6; nG <- 300L; muKa <- 0.45
  nm <- c("THETA[1]", "THETA[2]", "THETA[3]", "THETA[4]", "ETA[1]")

  e <- new.env()
  e$predOnly <- m
  e$rxControl <- rxode2::rxControl(atol=1e-10, rtol=1e-10)
  e$param <- stats::setNames(c(th, 0), nm)
  e$data  <- dat

  rxode2::rxSetSeed(11)
  est <- rpemEstepK1Draw(e, c(th, 0), 4L, matrix(om, 1, 1), nG, 1L)
  on.exit(rpemFree(), add=TRUE)

  etas <- as.numeric(est$eta)   # layout [(i-1)*nG + j]
  logp <- as.numeric(est$logp)

  # Deterministic IS reference: theta_ij = muKa + eta_ij; w_ij = softmax over j
  wList <- vector("list", nid); thetaList <- vector("list", nid); ESub <- numeric(nid)
  for (i in seq_len(nid)) {
    idx <- ((i - 1L) * nG + 1L):(i * nG)
    lp <- logp[idx]; w <- exp(lp - max(lp)); w <- w / sum(w)
    theta <- muKa + etas[idx]
    wList[[i]] <- w; thetaList[[i]] <- theta
    ESub[i] <- sum(theta * w)
  }
  muRef <- mean(ESub)
  omegaRef <- mean(vapply(seq_len(nid), function(i)
    sum(wList[[i]] * (thetaList[[i]] - muRef)^2), numeric(1)))

  rxode2::rxSetSeed(22)
  ms <- rpemMstepK1(muKa, 0.7, 400000L, 40000L)   # addSd0 = model add.sd

  expect_gt(ms$accept, 0); expect_lt(ms$accept, 1)
  expect_equal(ms$mu[1], muRef, tolerance = 0.02)
  expect_equal(ms$omega[1, 1], omegaRef, tolerance = 0.1)
})
