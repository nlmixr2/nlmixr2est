## Omega off-diagonal (correlated eta) estimation across the newer estimation
## methods.  A correlated 2-eta block `eta.cl + eta.v ~ c(...)` must be
## ESTIMATED, not silently frozen or dropped: the data are simulated with a
## STRONG positive eta.cl/eta.v correlation and each method must recover a
## substantively positive correlation (not merely "moved off the ini value",
## which a regression zeroing the off-diagonal would also satisfy).
## Real fits -> weekly slow batch.

nmTest({
  ## simulation truth: rho = 0.75 on omega diag 0.09
  .omTrue <- matrix(c(0.09, 0.0675, 0.0675, 0.09), 2, 2,
                    dimnames = list(c("eta.cl", "eta.v"), c("eta.cl", "eta.v")))

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
    f <- nlmixr2(.omCorMod, .omData, est = "advi",
                 control = adviControl(iters = 300L, print = 0L))
    .expectOffDiagEstimated(f)
    .idf <- f$iniDf
    .offRow <- .idf[!is.na(.idf$neta1) & .idf$neta1 != .idf$neta2, , drop = FALSE]
    expect_equal(nrow(.offRow), 1L)
    expect_equal(as.numeric(.offRow$est), f$omega[1L, 2L], tolerance = 1e-8)
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
