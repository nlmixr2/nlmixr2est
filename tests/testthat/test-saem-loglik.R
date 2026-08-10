nmTest({
  # General log-likelihood endpoints (ll() ~ expr) are supported by BOTH saem and
  # fsaem the saemix way: the model returns the per-observation log-likelihood as
  # its prediction, the observation loss is -ll, and the standard MCMC kernels run
  # unchanged (fsaem adds the IMH accelerator on top).

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

  # [pkpd.two.endpoint()] with the PD endpoint written as its own log-likelihood
  # instead of add(pd.sd), so the model mixes a normal and an ll() endpoint
  mixedNormLl <- function() {
    ini({
      tcl <- -3.2; tv <- -1; tec50 <- 0.5
      eta.cl ~ 0.09; eta.v ~ 0.09
      pk.sd <- 0.4; pd.sd <- 4
    })
    model({
      ka <- exp(0.5)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      e0 <- 100; ec50 <- exp(tec50)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      eff <- e0 * (1 - cp / (ec50 + cp))
      cp ~ add(pk.sd) | cp
      ll(eff) ~ -log(pd.sd) - 0.5*log(2*pi) - 0.5*((DV - eff)/pd.sd)^2
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

  test_that("fsaem fits the same general log-likelihood endpoint", {
    .d <- .mkTte(1L)
    .f <- suppressMessages(nlmixr2(expTte, .d, est = "fsaem",
      control = saemControl(nBurn = 150, nEm = 80, nmc = 3, seed = 1, print = 0L,
                            calcTables = FALSE)))
    expect_true(.f$saemControl$fast)
    expect_equal(exp(fixef(.f)[["tlam"]]), 40, tolerance = 0.2)
  })

  test_that("saem and fsaem agree on a general log-likelihood endpoint", {
    .d <- .mkTte(2L)
    .ctl <- saemControl(nBurn = 200, nEm = 100, nmc = 3, seed = 3, print = 0L,
                        calcTables = FALSE)
    .ss <- suppressMessages(nlmixr2(expTte, .d, est = "saem",  control = .ctl))
    .fs <- suppressMessages(nlmixr2(expTte, .d, est = "fsaem", control = .ctl))
    # both land at the same MLE (the observation likelihood is identical)
    expect_lt(abs(fixef(.ss)[["tlam"]] - fixef(.fs)[["tlam"]]), 0.08)
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

  test_that("the residual M-step skips an ll() endpoint of a mixed model (#873)", {
    # A mixed norm + ll() model carries the scalar distribution==1, so the
    # `if (distribution != 4)` guard on the residual M-step does not fire and the
    # loop used to run over the ll() endpoint as well.  Nothing accumulates into
    # statrese for that endpoint, so it stored sigma2 = 0 -- the value the FIM's
    # residual slot then divides by.
    .fitOn <- function(.d) {
      .f <- suppressMessages(nlmixr2(mixedNormLl, .d, est = "saem",
        control = saemControl(nBurn = 10, nEm = 5, nmc = 3, seed = 1, print = 0L,
                              calcTables = FALSE, covMethod = 0L)))
      .f$saem$res_info
    }
    .r1 <- .fitOn(mkPkpdTwoEndpointData(n = 12L))
    .r2 <- .fitOn(mkPkpdTwoEndpointData(n = 12L, seed = 7L))
    # endpoint 1 is the additive normal one, endpoint 2 the ll()
    expect_equal(as.integer(.r1$res_mod), c(1L, 0L))
    # the mechanism, without pinning saem's sigma2 initialization constant: the
    # ll() endpoint's slot is untouched by the M-step, so it cannot depend on the
    # data, while the normal endpoint's does.  The bug drove the ll() slot to
    # exactly 0 -- also data-independent, hence the second assertion.
    expect_equal(as.numeric(.r1$sigma2)[2], as.numeric(.r2$sigma2)[2])
    expect_gt(as.numeric(.r1$sigma2)[2], 0)
    expect_false(isTRUE(all.equal(as.numeric(.r1$sigma2)[1], as.numeric(.r2$sigma2)[1])))
    # and the normal endpoint is still estimated
    expect_true(is.finite(as.numeric(.r1$sigma2)[1]))
    expect_gt(as.numeric(.r1$sigma2)[1], 0)
  })
})
