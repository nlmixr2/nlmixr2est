nmTest({
  test_that("plateau criterion needs every parameter flat (#1395)", {
    .k <- 0:20
    .theta <- cbind(3.45 + 1e-5 * .k, 0.2 + 1e-5 * .k)
    .omFlat <- 0.1 + 1e-6 * .k
    .omDrift <- 0.1 * 1.01^.k
    .sig <- 0.7 + 1e-6 * .k
    expect_true(saemPlateauTest_(cbind(.theta, .omFlat, .sig), 2L, 1L, numeric(0), 1e-3))
    # theta flat but omega still drifting
    expect_false(saemPlateauTest_(cbind(.theta, .omDrift, .sig), 2L, 1L, numeric(0), 1e-3))
    # a theta near zero uses the absolute change; an omega never does
    .near0 <- 1e-4 * sin(.k)
    expect_true(saemPlateauTest_(cbind(.near0, .omFlat), 1L, 1L, numeric(0), 1e-3))
    expect_false(saemPlateauTest_(cbind(3.45, 0.001 + .near0), 1L, 1L, numeric(0), 1e-3))
    # an off-diagonal omega is scaled by sqrt(w_ii * w_jj), not by itself
    .off <- 1e-4 * cos(.k)
    expect_true(saemPlateauTest_(cbind(.theta, .omFlat, .off), 2L, 1L, 0.1, 1e-2))
    expect_false(saemPlateauTest_(cbind(.theta, .omFlat, .off), 2L, 1L, 0.001, 1e-2))
    # a non-finite value or a single row never plateaus
    expect_false(saemPlateauTest_(cbind(.theta, c(NA, .omFlat[-1])), 2L, 1L, numeric(0), 1e-3))
    expect_false(saemPlateauTest_(cbind(3.45, 0.1), 1L, 1L, numeric(0), 1e-3))
  })

  test_that("saemControl() validates the plateau options (#1395)", {
    .c <- saemControl(plateauStop = TRUE, plateauWindow = 20, plateauTol = 0.01, plateauMinFrac = 0.25)
    expect_true(.c$plateauStop)
    expect_identical(.c$plateauWindow, 20L)
    expect_equal(.c$plateauTol, 0.01)
    expect_equal(.c$plateauMinFrac, 0.25)
    expect_false(saemControl()$plateauStop)
    expect_error(saemControl(plateauStop = NA))
    expect_error(saemControl(plateauWindow = 0))
    expect_error(saemControl(plateauTol = -1))
    expect_error(saemControl(plateauMinFrac = 1.5))
  })

  test_that("plateau stop ends the EM phase like a smaller nEm (#1395)", {
    skip_on_cran()
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- log(2.7)
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .fit <- function(..., seed = 99) {
      suppressMessages(nlmixr2(
        one.cmt,
        nlmixr2data::theo_sd,
        "saem",
        control = saemControl(print = 0, nBurn = 100, nEm = 200, seed = seed, ...)
      ))
    }
    .plateauNote <- function(fit) {
      grep("plateaued", fit$runInfo, value = TRUE)
    }
    .stop <- function(fit) {
      as.integer(sub(".*iteration ([0-9]+) of.*", "\\1", .plateauNote(fit)))
    }
    .full <- .fit()
    expect_length(.plateauNote(.full), 0L)
    expect_equal(nrow(.full$parHistData), 300L)

    # a criterion that never holds leaves the fit untouched
    .never <- .fit(plateauStop = TRUE, plateauTol = 0)
    expect_length(.plateauNote(.never), 0L)
    expect_equal(.never$parFixedDf$Estimate, .full$parFixedDf$Estimate)
    expect_equal(.never$omega, .full$omega)
    expect_equal(.never$cov, .full$cov)

    .opt <- list(plateauStop = TRUE, plateauTol = 0.02, plateauWindow = 30L)
    .stopped <- do.call(.fit, .opt)
    .n <- .stop(.stopped)
    expect_length(.n, 1L)
    expect_true(.n >= 100L + 0.5 * 200L && .n < 300L)
    expect_lt(nchar(.plateauNote(.stopped)), 75L)
    # every per-iteration history ends at the stop
    expect_equal(nrow(.stopped$parHistData), .n)
    expect_equal(dim(.stopped$phiM)[3], .n)

    # the same fit as asking for nEm = .n - nBurn up front, covariance included
    .short <- suppressMessages(nlmixr2(
      one.cmt,
      nlmixr2data::theo_sd,
      "saem",
      control = saemControl(print = 0, nBurn = 100, nEm = .n - 100L)
    ))
    expect_equal(.stopped$parFixedDf$Estimate, .short$parFixedDf$Estimate)
    expect_equal(.stopped$omega, .short$omega)
    expect_equal(.stopped$cov, .short$cov)
    expect_equal(.stopped$objf, .short$objf)

    # deterministic under the seed and independent of the thread count
    .old <- rxode2::getRxThreads()
    on.exit(rxode2::setRxThreads(.old), add = TRUE)
    rxode2::setRxThreads(1L)
    .one <- do.call(.fit, .opt)
    expect_equal(.stop(.one), .n)
    expect_equal(.one$parFixedDf$Estimate, .stopped$parFixedDf$Estimate)
    rxode2::setRxThreads(.old)

    # early-stopped estimates agree with the full run within Monte Carlo error
    for (.seed in c(99, 2024)) {
      .f <- if (.seed == 99) .full else .fit(seed = .seed)
      .s <- if (.seed == 99) .stopped else do.call(.fit, c(.opt, list(seed = .seed)))
      expect_equal(.s$parFixedDf$Estimate, .f$parFixedDf$Estimate, tolerance = 0.05)
      expect_equal(diag(.s$omega), diag(.f$omega), tolerance = 0.1)
    }
  })
})
