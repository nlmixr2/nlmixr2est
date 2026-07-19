nmTest({
  # General log-likelihood endpoints (ll() ~ expr) are supported by saem the saemix
  # way: the model returns the per-observation log-likelihood as its prediction, the
  # observation loss is -ll, and the standard MCMC kernels run unchanged.

  # exponential time-to-event data with a subject random effect on the mean
  .mkTte <- function(seed = 1L, n = 150L, meanT = 40) {
    set.seed(seed)
    do.call(rbind, lapply(seq_len(n), function(i) {
      lami <- meanT * exp(rnorm(1, 0, sqrt(0.15)))
      data.frame(ID = i, TIME = lami * -log(runif(1)), DV = 1, EVID = 0, CMT = 1)
    }))
  }

  expTte <- function() {
    ini({ tlam <- log(25); eta.lam ~ 0.2 })
    model({
      lam <- exp(tlam + eta.lam)
      ll(dv) ~ -log(lam) - time / lam            # exponential event-time loglik
    })
  }

  test_that("saem fits a general log-likelihood (exponential TTE) endpoint", {
    .d <- .mkTte(1L)
    .f <- suppressMessages(nlmixr2(expTte, .d, est = "saem",
      control = saemControl(nBurn = 150, nEm = 80, nmc = 3, seed = 1, print = 0L,
                            calcTables = FALSE)))
    # recovers the population mean (true 40) from a poor start (25)
    expect_equal(exp(fixef(.f)[["tlam"]]), 40, tolerance = 0.2)
    # a random-effect variance was estimated
    expect_gt(.f$omega[1, 1], 0)
  })

  test_that("a general log-likelihood endpoint uses distribution=4 (no residual)", {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(expTte))
    # LL endpoint carries no residual bookkeeping
    expect_equal(.ui$saemResMod, 0L)
    expect_equal(.ui$saemModNumEst, 0L)
    # and the saem solve model emits the log-likelihood as rx_pred_ (rx_yj_ ~ 152)
    expect_true(any(grepl("rx_yj_", as.character(.ui$saemModel0))))
  })
})
