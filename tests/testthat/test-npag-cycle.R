## NPAG adaptive-grid cycle (npagCycle_) end-to-end on Theophylline: Sobol grid
## -> Psi -> Burke IPM -> condensation -> expansion -> convergence.  Runs a real
## (if modest) fit, so it lives in a weekly slow batch.

nmTest({
  test_that("npagCycle_ converges to a sane discrete distribution on theo_sd", {
    npMod <- function() {
      ini({ tka <- log(1.5); tv <- log(32); tke <- log(0.08)
        eta.ka ~ 0.3; eta.v ~ 0.1; eta.ke ~ 0.1
        add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); ke <- exp(tke + eta.ke)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - ke * center
        cp <- center / v
        cp ~ add(add.sd) })
    }
    ui <- rxode2::assertRxUi(npMod)
    ctl <- npagControl()
    dat <- nlmixr2data::theo_sd
    N <- length(unique(dat$ID))

    .npInnerSetup(ui, dat, matrix(0, N, 3L), ctl)
    on.exit(.npInnerFree(), add = TRUE)
    box <- .npEtaBox(ui, ctl)
    cores <- as.integer(rxode2::getRxThreads())

    r <- npagCycle_(box$lower, box$upper, points = 512L, cycles = 30L, cores = cores)

    # weights form a probability vector
    expect_equal(sum(r$weights), 1, tolerance = 1e-6)
    expect_true(all(r$weights >= -1e-10))
    # Caratheodory: the NPML has at most N support points
    nspp <- nrow(r$support)
    expect_true(nspp >= 1L && nspp <= N)
    expect_equal(ncol(r$support), 3L)
    # support points stay within the box
    for (j in 1:3) {
      expect_true(all(r$support[, j] >= box$lower[j] - 1e-8 &
                        r$support[, j] <= box$upper[j] + 1e-8))
    }
    # objective / -2LL are finite and in a sane range for this 12-subject PK model
    expect_true(is.finite(r$objf))
    m2ll <- -2 * r$objf
    expect_true(m2ll > 50 && m2ll < 5000)
    expect_true(r$converged || r$cycles == 30L)

    # the weighted-mean eta is a small deviation (etas are individual offsets)
    wmean <- colSums(r$support * r$weights)
    expect_true(all(abs(wmean) < 2))
  })
})
