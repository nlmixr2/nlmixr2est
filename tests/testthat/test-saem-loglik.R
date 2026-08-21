nmTest({
  # General log-likelihood endpoints (ll() ~ expr) are supported by saem the saemix
  # way: the model returns the per-observation log-likelihood as its prediction, the
  # observation loss is -ll, and the standard MCMC kernels run unchanged.

  # exponential time-to-event data with a subject random effect on the mean
  .mkTte <- function(seed = 1L, n = 150L, meanT = 40) {
    .testSeed(seed)
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

  test_that("a Gaussian twin agrees with the equivalent add() model (#871)", {
    # An ll() written as the exact normal log-density and the equivalent add()
    # model are the same likelihood, so they must give the same objective function
    # value and the same standard errors.  The residual SD is fixed in both so the
    # two fits estimate the identical free parameter set, and it is carried on the
    # log scale in the ll() because an SD written bare inside an ll() has no
    # positivity constraint.
    mAdd <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        add.sd <- fixed(0.7)
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    mLl <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        lsd <- fixed(log(0.7))
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        sd <- exp(lsd)
        cp <- linCmt()
        ll(err) ~ -lsd - 0.5 * log(2 * pi) - 0.5 * ((DV - cp) / sd)^2
      })
    }
    ctl <- saemControl(nBurn = 300, nEm = 400, seed = 42L, print = 0L,
                       covMethod = "linFim", calcTables = FALSE)
    fA <- .nlmixr(mAdd, theo_sd, est = "saem", control = ctl)
    fL <- .nlmixr(mLl,  theo_sd, est = "saem", control = ctl)

    # the ll() fit really did take the general-likelihood path
    expect_equal(fL$ui$saemResMod, 0L)
    expect_equal(fA$ui$saemResMod, 1L)
    # the two fits land on the same estimates, so any difference below is the
    # objective/covariance code rather than the fit
    expect_equal(unname(fixef(fL)[1:3]), unname(fixef(fA)[1:3]), tolerance = 0.01)
    expect_equal(unname(diag(fL$omega)), unname(diag(fA$omega)), tolerance = 0.05)

    # objective: equal, not merely equal up to a 0.5*log(2*pi) per observation.
    # Before the fix the ll() value was 766.5 against the add() model's 123.7.
    expect_equal(fL$objf, fA$objf, tolerance = 0.005)

    # standard errors: theta and the Omega variances.  Before the fix the ratios
    # ran from 4.0 to 80.5; the residual spread is the two fits' own difference.
    .k <- intersect(rownames(fL$cov), rownames(fA$cov))
    expect_true(all(c("tka", "tcl", "tv", "om.eta.ka", "om.eta.cl", "om.eta.v") %in% .k))
    expect_equal(unname(sqrt(diag(fL$cov))[.k]), unname(sqrt(diag(fA$cov))[.k]),
                 tolerance = 0.05)
  })

  test_that("a general log-likelihood endpoint uses distribution=4 (no residual)", {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(expTte))
    # LL endpoint carries no residual bookkeeping
    expect_equal(.ui$saemResMod, 0L)
    expect_equal(.ui$saemModNumEst, 0L)
    # and the saem solve model emits the log-likelihood as rx_pred_ (rx_yj_ ~ 152)
    expect_true(any(grepl("rx_yj_", as.character(.ui$saemModel0))))
  })

  # t()/cauchy() are not literal ll() syntax, but rxode2's own FOCEi line
  # generator reduces them to the same shape (rx_pred_ ~ an explicit
  # log-density, rx_r_ ~ 0) as ll() -- so saem should treat them the same way
  # rather than erroring ("t isn't supported yet" / "Distribution not
  # supported").  This does not assert convergence quality; that is covered
  # separately (see test-saem-loglik-focei.R).
  mT <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      add.sd <- 0.7
      nu <- fixed(8)
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      cp <- linCmt()
      cp ~ add(add.sd) + dt(nu)
    })
  }

  mCauchy <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      cp <- linCmt()
      cp ~ add(add.sd) + cauchy()
    })
  }

  test_that("t()/cauchy() endpoints dispatch as general-likelihood in saem (#Phase0)", {
    .uiT <- rxode2::rxUiDecompress(rxode2::rxode2(mT))
    expect_equal(.uiT$saemResMod, 0L)
    expect_equal(.uiT$saemModNumEst, 0L)
    .m0T <- as.character(.uiT$saemModel0)
    expect_true(any(grepl("llikT", .m0T)))
    expect_true(any(grepl("rx_r_ ~ 0", .m0T)))

    .uiC <- rxode2::rxUiDecompress(rxode2::rxode2(mCauchy))
    expect_equal(.uiC$saemResMod, 0L)
    expect_equal(.uiC$saemModNumEst, 0L)
    .m0C <- as.character(.uiC$saemModel0)
    expect_true(any(grepl("llikCauchy", .m0C)))
    expect_true(any(grepl("rx_r_ ~ 0", .m0C)))
  })

  test_that("t()/cauchy() endpoints fit without erroring in saem (#Phase0)", {
    ctl <- saemControl(nBurn = 20, nEm = 20, print = 0L, calcTables = FALSE)
    fT <- suppressWarnings(.nlmixr(mT, theo_sd, est = "saem", control = ctl))
    expect_equal(fT$ui$saemResMod, 0L)
    expect_true(is.finite(fT$objf))

    fC <- suppressWarnings(.nlmixr(mCauchy, theo_sd, est = "saem", control = ctl))
    expect_equal(fC$ui$saemResMod, 0L)
    expect_true(is.finite(fC$objf))
  })
})
