# Regression: a population parameter (theta) initialized at exactly 0 must still
# be estimated by FOCEi.
#
# The default nlmixr2 scaling constant is 1/|initPar|.  For a parameter
# initialized at 0 this is 1/0 = Inf, which clamps to scaleCmax and makes the
# parameter effectively unoptimizable -- a tiny step in the scaled space becomes
# an enormous step in the parameter, the line search rejects it, and the value
# stays frozen at ~0.  getScaleC() now falls back to unit scaling when initPar
# is 0 so the parameter is estimated normally.

nmTest({
  test_that("a theta initialized at exactly 0 is estimated, not frozen", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")

    # simulate a 1-cmt oral model with a strong linear body-weight effect on CL
    .sim <- function() {
      ini({
        tka <- 0.45; tcl <- -1.0; tv <- 3.45; b.cl <- 0.03
        eta.ka ~ 0.3; eta.cl ~ 0.02; eta.v ~ 0.05; add.sd <- 0.15
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + b.cl * WT + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    set.seed(11)
    .obsT <- seq(0.5, 60, length.out = 12)
    .ev <- do.call(rbind, lapply(seq_len(60), function(id) {
      wt <- round(runif(1, 50, 110))
      e <- as.data.frame(rxode2::et(amt = 320) |> rxode2::et(.obsT))
      e$ID <- id
      e$WT <- wt
      e
    }))
    .s <- rxode2::rxSolve(.sim(), .ev, returnType = "data.frame")
    .dose <- do.call(rbind, lapply(seq_len(60), function(id) {
      data.frame(ID = id, time = 0, DV = NA, amt = 320, evid = 1,
                 WT = .ev$WT[.ev$ID == id][1])
    }))
    .obs <- data.frame(ID = .s$id, time = .s$time, DV = .s$sim, amt = NA,
                       evid = 0, WT = .s$WT)
    .d <- rbind(.dose, .obs)
    .d <- .d[order(.d$ID, .d$time, -.d$evid), ]

    # candidate model: covariate coefficient `b` initialized at exactly 0
    .cand <- function() {
      ini({
        tka <- 0.45; tcl <- 1.0; tv <- 3.45; b <- 0
        eta.ka ~ 0.3; eta.cl ~ 0.1; eta.v ~ 0.1; add.sd <- 0.3
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + b * WT)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    fit <- suppressWarnings(
      nlmixr2(.cand, .d, "focei", foceiControl(print = 0L, covMethod = ""))
    )

    b <- unname(fit$parFixedDf["b", "Estimate"])
    # before the fix `b` stayed pinned at ~ -2e-4 (the finite-difference step);
    # it should now move well away from 0 toward the true value (0.03)
    expect_true(abs(b) > 0.01)
    # and land in a sensible neighborhood of the simulated effect (0.03)
    expect_true(b > 0.015 && b < 0.05)
  })
})
