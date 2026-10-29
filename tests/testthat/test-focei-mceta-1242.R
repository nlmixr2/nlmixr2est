nmTest({
  # #1242: mceta >= 1 ranked its Omega draws by the inner objective at the
  # STARTING point and optimized only the winner (plus eta=0), so a draw that
  # would have converged to a lower optimum was never optimized.  Every draw is
  # now optimized and the best converged result kept.
  #
  # A flip-flop one-compartment model has two optima per subject (ka and ke
  # swap, with V rescaled).  At these parameters eta=0 converges to the higher
  # one, while the lower one lies in the basin of many Omega draws.
  .ffMod <- function() {
    ini({
      tka <- log(0.4)
      tke <- log(0.5)
      tv <- log(6)
      eta.ka ~ 0.5
      eta.ke ~ 0.5
      eta.v ~ 0.5
      prop.sd <- 0.05
    })
    model({
      ka <- exp(tka + eta.ka)
      v <- exp(tv + eta.v)
      cl <- exp(tke + eta.ke) * v
      linCmt() ~ prop(prop.sd)
    })
  }

  .ffData <- function() {
    .t <- c(0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24)
    .c <- 100 / (10 * (1 - 0.2)) * (exp(-0.2 * .t) - exp(-.t))
    do.call(rbind, lapply(1:12, function(id) {
      rbind(
        data.frame(ID = id, TIME = 0, DV = NA, EVID = 1, AMT = 100),
        data.frame(ID = id, TIME = .t, DV = .c, EVID = 0, AMT = 0)
      )
    }))
  }

  .ffFit <- function(dat, mceta, etaMat = NULL) {
    suppressWarnings(suppressMessages(
      nlmixr2(
        .ffMod,
        dat,
        "focei",
        foceiControl(
          print = 0L,
          covMethod = "",
          maxOuterIterations = 0L,
          calcTables = FALSE,
          mceta = mceta,
          etaMat = etaMat
        )
      )
    ))
  }

  test_that("mceta optimizes every draw and keeps the lowest converged optimum", {
    skip_on_cran()
    .dat <- .ffData()
    # Every subject started inside the lower basin.
    .low <- matrix(c(0.91, -0.915, 0.51), 12, 3, byrow = TRUE)
    .fLow <- .ffFit(.dat, -2L, etaMat = .low)
    expect_true(all(.fLow$eta$eta.ka > 0))

    # From eta=0 every subject converges to the higher optimum.
    .f0 <- .ffFit(.dat, 0L)
    expect_true(all(.f0$eta$eta.ka < 0))
    expect_gt(.f0$objf, .fLow$objf + 1)

    # Ranking the draws at their starting point left subjects in the higher
    # optimum; optimizing each draw finds the lower one for all of them.
    .f10 <- .ffFit(.dat, 10L)
    expect_true(all(.f10$eta$eta.ka > 0))
    expect_equal(.f10$objf, .fLow$objf, tolerance = 1e-4)
    expect_gt(.f10$env$nMcetaStart[["sample"]], 0L)
  })
})
