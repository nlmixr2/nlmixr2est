# End-to-end RPEM K=1 (E-step -> conjugate M-step -> repeat) must agree with FOCEI
# on the same well-identified data (design/rpem/11 Bar 2, criteria C2.2).  Fixed
# effects (tcl, tv, add.sd) are held at truth in both fits since RPEM's numeric
# fixed-effect update is not wired yet; RPEM estimates mu (typical log-ka) and Omega.

test_that("RPEM E-M loop agrees with FOCEI on population mu and Omega (K=1)", {
  skip_on_cran()
  skip_on_ci()  # heavy: a FOCEI fit plus a multi-iteration RPEM loop

  struct <- rxode2::rxode2({ ka <- exp(tka+eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(7)
  nsub <- 40L; trueTka <- 0.45; trueOm <- 0.3; addSd <- 0.15; obsT <- seq(0.25, 24, by=1.5)
  etasTrue <- rnorm(nsub, 0, sqrt(trueOm))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt=100, cmt="depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params=c(tka=trueTka, tcl=1, tv=3.45, eta=etasTrue[i]),
                         events=ev, returnType="data.frame", addDosing=FALSE)
    d <- as.data.frame(ev); d$id <- i; d$DV <- 0; o <- d$evid == 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, addSd); d
  }))

  # FOCEI reference: fix tcl, tv at truth (RPEM holds those); estimate tka, Omega,
  # and add.sd (RPEM now estimates the additive residual too).
  mod <- function() {
    ini({ tka <- 0.1; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.3; eta.ka ~ 1.0 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  fit <- suppressMessages(nlmixr2(mod, dat, est="focei", control=foceiControl(print=0)))
  fTka <- fit$parFixedDf["tka", "Estimate"]; fOm <- fit$omega[1, 1]
  fSd <- fit$parFixedDf["add.sd", "Estimate"]

  # RPEM E-M loop.
  one.cmt <- function() {
    ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; add.sd <- 0.15; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  m <- rxode2::rxode2(one.cmt)$rpemRxModel$predOnly
  nm <- c("THETA[1]", "THETA[2]", "THETA[3]", "THETA[4]", "ETA[1]")
  e <- new.env(); e$predOnly <- m; e$rxControl <- rxode2::rxControl(atol=1e-8, rtol=1e-8)
  e$param <- stats::setNames(c(trueTka, 1, 3.45, addSd, 0), nm); e$data <- dat

  mu <- 0.1; Om <- 1.0; sdHat <- 0.3; nG <- 300L; niter <- 30L; collect <- 12L
  base <- c(mu, 1, 3.45, sdHat, 0)
  muTr <- numeric(niter); omTr <- numeric(niter); sdTr <- numeric(niter)
  for (it in seq_len(niter)) {
    base[1] <- mu; base[4] <- sdHat            # estimate mu, Omega, AND add.sd
    rxode2::rxSetSeed(2000 + it)
    rpemEstepK1Draw(e, base, 4L, matrix(Om, 1, 1), nG, 1L)
    ms <- rpemMstepK1(mu, sdHat, 80000L, 8000L)
    mu <- ms$mu[1]; Om <- ms$omega[1, 1]; sdHat <- ms$addSd
    muTr[it] <- mu; omTr[it] <- Om; sdTr[it] <- sdHat
  }
  rpemFree()

  # RPEM estimate = mean over the last `collect` (converged) iterations.
  muHat <- mean(tail(muTr, collect)); omHat <- mean(tail(omTr, collect))
  sdEst <- mean(tail(sdTr, collect))
  expect_equal(muHat, fTka, tolerance = 0.04)
  expect_equal(omHat, fOm,  tolerance = 0.10)
  expect_equal(sdEst, fSd,  tolerance = 0.10)
})
