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

  test_that("saemControl(phi1Hessian=TRUE) fits end-to-end on an ODE model (#Phase4)", {
    # An ODE (non-linCmt()) general-likelihood endpoint takes the ANALYTIC
    # Hessian path (odeSlotHess2, which carries its own event-sensitivity
    # shape distinct from odeSlotPred's "none") -- this is the only test that
    # exercises phi1Hessian=TRUE through an actual fit rather than a
    # standalone model build (test-saem-phi1-inner.R), which is what let a
    # real bug (odeSlotHess2 solved under phi1Objective's OpenMP region with
    # no OdeSwapEsBatch installed outside it, corrupting handle_evid's
    # scratch) ship undetected.
    mLl2 <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        lsd <- fixed(log(0.7))
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        sd <- exp(lsd)
        ll(err) ~ -lsd - 0.5 * log(2 * pi) - 0.5 * ((DV - cp) / sd)^2
      })
    }
    .n0 <- saemPhi1RefineN_()
    ctl <- saemControl(nBurn = 40, nEm = 40, nmc = 3, seed = 42L, print = 0L,
                       covMethod = "", calcTables = FALSE, phi1Hessian = TRUE)
    f <- suppressWarnings(.nlmixr(mLl2, theo_sd, est = "saem", control = ctl))
    expect_true(is.finite(f$objf))
    expect_equal(unname(fixef(f)[["tka"]]), 0.45, tolerance = 0.2)
    # confirms the analytic-Hessian path actually ran, not just that the fit
    # returned something
    expect_gt(saemPhi1RefineN_(), .n0)
  })

  test_that("a covariate-on-a-mu-ref-theta Gaussian twin still agrees with add() (#Phase6)", {
    # .saemPhi1TargetMap declines a covariate on a mu-referenced theta (v1
    # scope, see test-saem-phi1-inner.R), so this twin's phi1 step falls all
    # the way back to the historic Robbins-Monro SA-recursion for every
    # theta -- exercising that the decline-and-fallback path (not Phase 4's
    # new step) still lands on the right answer for a general-likelihood
    # endpoint, broadening the #871 twin beyond the etas-only case.
    mAdd <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45; tka.wt <- 0.01
        add.sd <- fixed(0.7)
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + tka.wt * WT + eta.ka)
        cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    mLl <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45; tka.wt <- 0.01
        lsd <- fixed(log(0.7))
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + tka.wt * WT + eta.ka)
        cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        sd <- exp(lsd)
        cp <- linCmt()
        ll(err) ~ -lsd - 0.5 * log(2 * pi) - 0.5 * ((DV - cp) / sd)^2
      })
    }
    ctl <- saemControl(nBurn = 300, nEm = 400, seed = 42L, print = 0L,
                       covMethod = "", calcTables = FALSE)
    fA <- .nlmixr(mAdd, theo_sd, est = "saem", control = ctl)
    fL <- .nlmixr(mLl,  theo_sd, est = "saem", control = ctl)

    # this twin really does decline the new phi1 step (confirms the premise
    # of this test, not just its conclusion)
    expect_false(isTRUE(fL$ui$saemPhi1Inner$ok))

    expect_equal(unname(fixef(fL)[c("tka", "tcl", "tv", "tka.wt")]),
                 unname(fixef(fA)[c("tka", "tcl", "tv", "tka.wt")]),
                 tolerance = 0.05)
    expect_equal(unname(diag(fL$omega)), unname(diag(fA$omega)), tolerance = 0.1)
    expect_equal(fL$objf, fA$objf, tolerance = 0.02)
  })

  # Shared 2-endpoint data for the multi-endpoint general-likelihood twins
  # below: "cp" and a second, independently-noised observation "cp2" of the
  # same underlying concentration (deterministically offset so it is a
  # distinct condition, not literally the same numbers).
  .mkTwoEp <- function(seed = 1L, nSub = 20L) {
    .testSeed(seed)
    simMod <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
            eta.ka ~ 0.1; eta.cl ~ 0.1; eta.v ~ 0.05 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
              cp <- center / v
              cp ~ add(0.3) })
    }
    simUi <- rxode2::rxUiDecompress(rxode2::rxode2(simMod))
    ev <- rxode2::et(amt = 320, time = 0) |> rxode2::et(seq(0.5, 24, by = 2))
    sim <- rxode2::rxSolve(simUi, ev, params = c(tka = 0.45, tcl = 1, tv = 3.45), nSub = nSub,
                            omega = lotri::lotri(eta.ka ~ 0.1, eta.cl ~ 0.1, eta.v ~ 0.05))
    df <- as.data.frame(sim); names(df)[names(df) == "sim.id"] <- "id"
    d1 <- data.frame(id = df$id, time = df$time, amt = NA_real_, evid = 0,
                      dv = pmax(df$cp + rnorm(nrow(df), 0, 0.3), 0.01), cmt = "cp")
    d2 <- data.frame(id = df$id, time = df$time, amt = NA_real_, evid = 0,
                      dv = pmax(df$cp * 1.3 + rnorm(nrow(df), 0, 0.4), 0.01), cmt = "cp2")
    dose <- data.frame(id = sort(unique(df$id)), time = 0, amt = 320, evid = 101,
                        dv = NA_real_, cmt = "depot")
    dat <- rbind(dose, d1, d2); dat[order(dat$id, dat$time), ]
  }
  .twoEpCtl <- saemControl(nBurn = 100, nEm = 100, nmc = 3, seed = 42, print = 0L,
                           covMethod = "", calcTables = FALSE)
  .twoEpGauss <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
          add.pk1 <- 0.3; add.pk2 <- 0.4
          eta.ka ~ 0.1; eta.cl ~ 0.1; eta.v ~ 0.05 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
            d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
            cp <- center / v; cp2 <- cp * 1.3
            cp  ~ add(add.pk1) | cp
            cp2 ~ add(add.pk2) | cp2 })
  }

  test_that("a multi-endpoint Gaussian twin agrees with all-ll() endpoints", {
    # .saemGeneralLik() used to require exactly one prediction endpoint
    # (length(predDf$cond) == 1L), so a model with SEVERAL general-likelihood
    # endpoints (every condition non-"norm") was silently misclassified as
    # distribution="normal" instead of "general" -- do_mcmc's distribution==4
    # branch treats each row's rx_pred_ as -loglik directly regardless of
    # which endpoint condition produced it, so this needed no per-endpoint
    # C++ change, only widening the R-side gate to length(cond) >= 1L with
    # any(distribution != "norm").
    dat <- .mkTwoEp(1L)
    pkLl <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
            lsd1 <- log(0.3); lsd2 <- log(0.4)
            eta.ka ~ 0.1; eta.cl ~ 0.1; eta.v ~ 0.05 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
              cp <- center / v; cp2 <- cp * 1.3
              sd1 <- exp(lsd1); sd2 <- exp(lsd2)
              ll(cp)  ~ -0.5 * log(2 * pi) - lsd1 - 0.5 * ((DV - cp) / sd1)^2
              ll(cp2) ~ -0.5 * log(2 * pi) - lsd2 - 0.5 * ((DV - cp2) / sd2)^2 })
    }
    fA <- .nlmixr(.twoEpGauss, dat, est = "saem", control = .twoEpCtl)
    fL <- .nlmixr(pkLl,        dat, est = "saem", control = .twoEpCtl)

    # the ll() twin really did take the general (distribution=4) path, over
    # every endpoint -- confirms the premise, not just the conclusion
    expect_true(.saemGeneralLik(fL$ui))
    expect_equal(length(fL$ui$predDf$cond), 2L)

    expect_equal(unname(fixef(fL)[c("tka", "tcl", "tv")]),
                 unname(fixef(fA)[c("tka", "tcl", "tv")]),
                 tolerance = 0.02)
  })

  test_that("a mixed norm+ll() multi-endpoint twin agrees with its Gaussian equivalent", {
    # A "norm" condition mixed with a genuine general-likelihood condition
    # (some endpoints "norm", some not) used to fall all the way back to
    # distribution="normal" for the WHOLE fit, silently mis-scoring the
    # general-likelihood endpoint's rx_pred_ (a log-density) as if it were a
    # plain prediction.  Fixed by promoting every "norm" condition to
    # "dnorm" for dispatch purposes whenever it is mixed with a genuine
    # general-likelihood condition (.createSaemLineObject,
    # R/saemRxUiGetModel.R) -- llikNorm() is the exact same normal
    # log-density, just expressed as a proper general likelihood, mirroring
    # FOCEi's own rxUiGet.predDfFocei() promotion (R/focei.R). This needed no
    # C++ change: do_mcmc's distribution==4 branch was already
    # endpoint-agnostic.
    dat <- .mkTwoEp(1L)
    pkMix <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
            add.pk1 <- 0.3; lsd2 <- log(0.4)
            eta.ka ~ 0.1; eta.cl ~ 0.1; eta.v ~ 0.05 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
              cp <- center / v; cp2 <- cp * 1.3
              sd2 <- exp(lsd2)
              cp ~ add(add.pk1) | cp
              ll(cp2) ~ -0.5 * log(2 * pi) - lsd2 - 0.5 * ((DV - cp2) / sd2)^2 })
    }
    fA <- .nlmixr(.twoEpGauss, dat, est = "saem", control = .twoEpCtl)
    fM <- .nlmixr(pkMix,       dat, est = "saem", control = .twoEpCtl)

    # the mixed fit really did take the general (distribution=4) path, and
    # its "norm" condition's residual bookkeeping was correctly suppressed
    # (promoted to dnorm, not left on the closed-form M-step)
    expect_true(.saemGeneralLik(fM$ui))
    expect_equal(unname(fM$ui$saemResMod), c(0L, 0L))

    expect_equal(unname(fixef(fM)[c("tka", "tcl", "tv")]),
                 unname(fixef(fA)[c("tka", "tcl", "tv")]),
                 tolerance = 0.02)
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
    # regression guard: t()/cauchy()'s residual-family theta (add.sd, nu) used
    # to be completely undeclared for a general-likelihood endpoint (absent
    # from params()/init/fixed bookkeeping), and separately, a FIXED one was
    # reported back from a stale residual-M-step slot (.resMat) instead of the
    # kernel's own tracked value -- a fixed theta must come back EXACTLY fixed.
    expect_equal(unname(fixef(fT)[["nu"]]), 8)

    fC <- suppressWarnings(.nlmixr(mCauchy, theo_sd, est = "saem", control = ctl))
    expect_equal(fC$ui$saemResMod, 0L)
    expect_true(is.finite(fC$objf))
  })

  test_that("a fixed general-likelihood residual-family theta reports its true kernel value (#Phase0)", {
    # dt(nu)'s add.sd is also fixed here (isolating this test from any
    # phi1/eta convergence noise) -- fixef() must match the kernel's own
    # Plambda, not a stale value from the (unused, res.mod==0) residual M-step.
    mTFixed <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        add.sd <- fixed(0.7)
        nu <- fixed(8)
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        cp <- linCmt()
        cp ~ add(add.sd) + dt(nu)
      })
    }
    ctl <- saemControl(nBurn = 20, nEm = 20, seed = 1L, print = 0L, calcTables = FALSE)
    fT <- suppressWarnings(.nlmixr(mTFixed, theo_sd, est = "saem", control = ctl))
    expect_equal(unname(fixef(fT)[["add.sd"]]), 0.7)
    expect_equal(unname(fixef(fT)[["nu"]]), 8)
    expect_equal(unname(fixef(fT)[["add.sd"]]), unname(fT$saem$Plambda[4, 1]))
  })
})
