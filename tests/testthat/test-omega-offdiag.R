## Omega off-diagonal (correlated eta) estimation across the newer estimation
## methods.  A correlated 2-eta block `eta.cl + eta.v ~ c(...)` must be
## ESTIMATED, not silently frozen or dropped: the data are simulated with a
## STRONG positive eta.cl/eta.v correlation and each method must recover a
## substantively positive correlation (not merely "moved off the ini value",
## which a regression zeroing the off-diagonal would also satisfy).
## Real fits -> weekly slow batch.

nmTest({
  ## simulation truth: rho = 0.75 (0.0675 on a 0.09 variance)
  .omSimMod <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.cl + eta.v ~ c(0.09,
                         0.0675, 0.09)
      add.sd <- 0.3 })
    model({ ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) })
  }

  ## estimation model: same structure, started at a weak correlation so a
  ## frozen/dropped off-diagonal cannot look like a recovered one
  .omCorMod <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.cl + eta.v ~ c(0.1,
                         0.01, 0.1)
      add.sd <- 0.7 })
    model({ ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) })
  }

  ## 60 subjects on the theo_sd sampling schedule
  .omData <- local({
    .ev <- nlmixr2data::theo_sd[nlmixr2data::theo_sd$ID == 1,
                                c("TIME", "AMT", "EVID", "DV")]
    .et <- rxode2::et(.ev)
    .s <- rxode2::rxSolve(.omSimMod, .et, nSub = 60L, seed = 1042L,
                          addDosing = TRUE, returnType = "data.frame")
    .d <- data.frame(ID = .s$sim.id, TIME = .s$time, AMT = .s$amt,
                     EVID = .s$evid, DV = .s$sim)
    .d$AMT[is.na(.d$AMT)] <- 0
    .d$DV[.d$EVID != 0] <- NA_real_
    ## a dose row and an observation share time 0; keep only the dose there
    .d <- .d[!(.d$TIME == 0 & .d$EVID == 0), , drop = FALSE]
    .d
  })

  ## a method that estimates the block must recover a clearly positive
  ## correlation; one that freezes (0.01/sqrt(.1*.1) = 0.1) or drops it (0)
  ## cannot clear this bar
  .expectOffDiagEstimated <- function(f) {
    .om <- f$omega
    expect_equal(dim(.om), c(2L, 2L))
    expect_equal(.om[1L, 2L], .om[2L, 1L])
    expect_true(all(is.finite(.om)))
    expect_true(all(eigen(.om, symmetric = TRUE, only.values = TRUE)$values > 0))
    .rho <- .om[1L, 2L] / sqrt(.om[1L, 1L] * .om[2L, 2L])
    expect_gt(.rho, 0.4)
    expect_lt(.rho, 1)
  }

  test_that("est='vae' estimates the omega off-diagonal", {
    f <- nlmixr2(.omCorMod, .omData, est = "vae",
                 control = vaeControl(itersBurnIn = 50L, iters = 100L,
                                      covariateSelection = FALSE, print = 0L))
    .expectOffDiagEstimated(f)
    ## the updated model carries the whole block: the off-diagonal iniDf row
    ## holds the estimate, symmetric with $omega
    .idf <- f$iniDf
    .offRow <- .idf[!is.na(.idf$neta1) & .idf$neta1 != .idf$neta2, , drop = FALSE]
    expect_equal(nrow(.offRow), 1L)
    expect_equal(as.numeric(.offRow$est), f$omega[1L, 2L], tolerance = 1e-8)
  })

  test_that("est='advi' estimates the omega off-diagonal", {
    f <- nlmixr2(.omCorMod, .omData, est = "emvi",
                 control = emviControl(iters = 300L, print = 0L))
    .expectOffDiagEstimated(f)
    .idf <- f$iniDf
    .offRow <- .idf[!is.na(.idf$neta1) & .idf$neta1 != .idf$neta2, , drop = FALSE]
    expect_equal(nrow(.offRow), 1L)
    expect_equal(as.numeric(.offRow$est), f$omega[1L, 2L], tolerance = 1e-8)
  })

  ## a FIXED block must come back exactly at its ini values: the M-step holds
  ## fixed entries, so a regression that estimates them anyway is caught here
  .omFixedMod <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.cl + eta.v ~ fixed(0.1,
                             0.01, 0.1)
      add.sd <- 0.7 })
    model({ ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) })
  }

  test_that("perNoCor holds the correlation, and releasing it estimates one", {
    ## The hold must ENGAGE (perNoCor = 1 -> the off-diagonal never moves off
    ## ini) and RELEASE (perNoCor = 0 -> it is estimated from the first M-step).
    ## Without both halves a regression that never lifts the hold would look
    ## exactly like "this model has no correlation".
    ## gammaIter >= iters ON PURPOSE: the hold is
    ## round(perNoCor * min(gammaIter, iters)) iterations, so with a SMALLER
    ## gammaIter even perNoCor = 1 leaves a tail that estimates the correlation
    ## (gammaIter = 40 of 60 iterations left 20 free and returned 0.0485, not 0).
    .ctl <- function(p) vaeControl(itersBurnIn = 20L, iters = 60L, gammaIter = 60L,
                                   perNoCor = p, covariateSelection = FALSE,
                                   print = 0L)
    ## perNoCor = 1 holds throughout: a FREE correlation is held at ZERO (saem's
    ## diagmat() rule), NOT at its ini value.  Holding it at ini while the
    ## variances move is what left the block non-positive-definite.
    fHold <- nlmixr2(.omCorMod, .omData, est = "vae", control = .ctl(1))
    expect_equal(unname(fHold$omega[1L, 2L]), 0, tolerance = 1e-10)
    fFree <- nlmixr2(.omCorMod, .omData, est = "vae", control = .ctl(0))
    expect_gt(abs(fFree$omega[1L, 2L] - 0.01), 1e-6)
    ## and the default hold still leaves room to estimate on a short run
    ## (nbCorrel is a fraction of min(gammaIter, iters), not of gammaIter)
    fDef <- nlmixr2(.omCorMod, .omData, est = "vae",
                    control = vaeControl(itersBurnIn = 20L, iters = 60L,
                                         covariateSelection = FALSE, print = 0L))
    expect_gt(abs(fDef$omega[1L, 2L] - 0.01), 1e-6)
  })

  test_that("a FIXED omega block is held by vae and advi", {
    fV <- nlmixr2(.omFixedMod, .omData, est = "vae",
                  control = vaeControl(itersBurnIn = 20L, iters = 40L,
                                       covariateSelection = FALSE, print = 0L))
    expect_equal(unname(fV$omega[1L, 2L]), 0.01, tolerance = 1e-10)
    expect_equal(unname(diag(fV$omega)), c(0.1, 0.1), tolerance = 1e-10)
    fA <- nlmixr2(.omFixedMod, .omData, est = "emvi",
                  control = emviControl(iters = 50L, print = 0L))
    expect_equal(unname(fA$omega[1L, 2L]), 0.01, tolerance = 1e-10)
    expect_equal(unname(diag(fA$omega)), c(0.1, 0.1), tolerance = 1e-10)
  })

  test_that("est='npag' estimates the omega off-diagonal", {
    f <- nlmixr2(.omCorMod, .omData, est = "npag",
                 control = npagControl(points = 256L, cycles = 15L,
                                       gammaOptimize = FALSE))
    .expectOffDiagEstimated(f)
  })

  test_that("est='npb' estimates the omega off-diagonal", {
    f <- nlmixr2(.omCorMod, .omData, est = "npb",
                 control = npbControl(points = 50L, burnin = 100L, nsamp = 100L,
                                      seed = 42L))
    .expectOffDiagEstimated(f)
  })
})
