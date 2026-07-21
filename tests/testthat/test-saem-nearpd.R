nmTest({
  # A dense BLOCK(4) BSV spread across a PK endpoint and two covariate endpoints
  # drives Omega (Gamma2_phi1) to a non-positive-definite matrix mid-run.  SAEM
  # used to abort with "inv_sympd(): matrix is singular or not positive
  # definite"; it now projects Omega to the nearest positive-definite matrix
  # (nmNearPD), warns once, and finishes.

  .mkDat <- function() {
    .testSeed(42)
    nid <- 30
    cov1 <- rnorm(nid, 70, 10)
    cov2 <- 0.5 * cov1 + rnorm(nid, 0, 5)
    pk <- do.call(rbind, lapply(seq_len(nid), function(i) data.frame(
      ID = i, TIME = c(0, 0.5, 1, 2, 4, 8), EVID = c(1, 0, 0, 0, 0, 0),
      AMT = c(100, 0, 0, 0, 0, 0), CMT = 1,
      DV = c(0, 5, 7, 6, 4, 2) * exp(rnorm(1, 0, 0.2)) + rnorm(6, 0, 0.2))))
    pk$CMT[pk$EVID == 0] <- 2
    c1 <- data.frame(ID = seq_len(nid), TIME = 0, EVID = 0, AMT = 0, CMT = 3, DV = cov1)
    c2 <- data.frame(ID = seq_len(nid), TIME = 0, EVID = 0, AMT = 0, CMT = 4, DV = cov2)
    dat <- rbind(pk, c1, c2)
    dat[order(dat$ID, dat$TIME, dat$CMT), ]
  }

  mod <- function() {
    ini({
      tcl <- 1; tv <- 3.5; tc1 <- 70; tc2 <- 40
      eta.cl + eta.v + eta.c1 + eta.c2 ~
        c(0.1,
          0.03, 0.1,
          0.03, 0.03, 100,
          0.02, 0.02, 0.5, 30)          # dense BLOCK(4) across endpoints
      add.sd <- 0.5
      c.sd1 <- fixed(0.01); c.sd2 <- fixed(0.01)
    })
    model({
      cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
      c1p <- tc1 + eta.c1;  c1p ~ add(c.sd1) | cov1
      c2p <- tc2 + eta.c2;  c2p ~ add(c.sd2) | cov2
    })
  }

  test_that("saem recovers from a non-positive-definite Omega via nearPD", {
    .d <- .mkDat()
    .f <- suppressWarnings(suppressMessages(
      nlmixr2(mod, .d, est = "saem",
              control = saemControl(nBurn = 50, nEm = 80, seed = 42, print = 0L))))
    # the fit completes instead of erroring ...
    expect_true(inherits(.f, "nlmixr2FitData"))
    # ... and the nearPD projection is actually exercised (mechanism-used check):
    # the one-time warning is recorded in the fit's runInfo.
    expect_true(any(grepl("positive definite", .f$runInfo, ignore.case = TRUE)))
  })
})
