# End-to-end fits exercising the analytic FOCEi machinery for general log-likelihood
# (ll()) and generalized (Poisson) endpoints: the exact-Hessian objective
# (H = Omega^-1 - sum d2(logLik)/deta2) and the batched-fd2 Almquist outer gradient.
# A fast=TRUE fit converges to the same MLE as the finite-difference fit.  Slow (multiple
# fits), so weekly-batched via .slowBatches in tests/testthat.R -- do NOT add skip_on_ci().

nmTest({

  .ll_ode <- function() {
    ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
          eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
            d/dt(depot)  <- -ka * depot
            d/dt(center) <-  ka * depot - cl / v * center
            cp <- center / v
            ll(err) ~ -0.5 * log(2 * pi) - log(add.sd) - 0.5 * ((DV - cp) / add.sd)^2 })
  }
  .ll_lincmt <- function() {
    ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
          eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
            cp <- linCmt()
            ll(err) ~ -0.5 * log(2 * pi) - log(add.sd) - 0.5 * ((DV - cp) / add.sd)^2 })
  }

  test_that("ll() fast (analytic) fit matches the finite-difference fit", {
    skip_on_cran(); skip_if_not_installed("nlmixr2data")
    d <- nlmixr2data::theo_sd
    ## cached fast=FALSE reference -- see helper-gradref.R
    fF <- .numRef("fit-fd-ll-ode", function() {
      .f <- suppressMessages(nlmixr2(.ll_ode, d, "focei",
              foceiControl(print = 0L, covMethod = "", fast = FALSE)))
      list(objf = .f$objf, theta = unname(.f$theta))
    })
    fT <- suppressMessages(nlmixr2(.ll_ode, d, "focei", foceiControl(print = 0L, covMethod = "", fast = TRUE)))
    # fast=TRUE for an ll() endpoint uses the exact-Hessian objective + analytic outer
    # gradient; it converges to the same MLE within optimizer tolerance
    expect_equal(as.numeric(fT$objf), as.numeric(fF$objf), tolerance = 1e-2)
    expect_equal(unname(fT$theta), fF$theta, tolerance = 1e-2)
    # ...and it got there on the all-C++ gradient.  Without this the fit could quietly
    # decline every evaluation to finite differences and still pass the two above.
    expect_gt(fT$env$nAnalyticGradDirect, 0L)
  })

  test_that("ll() with linCmt() falls back gracefully to finite differences", {
    skip_on_cran(); skip_if_not_installed("nlmixr2data")
    # linCmt() has no 2nd-order state sensitivities (solved form), so both the exact-Hessian
    # objective and the analytic outer gradient are out of scope; fast=TRUE must transparently
    # fall back to finite differences and still converge to the same MLE (not error).
    d <- nlmixr2data::theo_sd
    ## cached fast=FALSE reference -- see helper-gradref.R
    fF <- .numRef("fit-fd-ll-lincmt", function() {
      .f <- suppressMessages(nlmixr2(.ll_lincmt, d, "focei",
              foceiControl(print = 0L, covMethod = "", fast = FALSE)))
      list(objf = .f$objf, theta = unname(.f$theta))
    })
    fT <- suppressMessages(nlmixr2(.ll_lincmt, d, "focei", foceiControl(print = 0L, covMethod = "", fast = TRUE)))
    expect_equal(as.numeric(fT$objf), as.numeric(fF$objf), tolerance = 1e-2)
    expect_equal(unname(fT$theta), fF$theta, tolerance = 1e-2)
    expect_equal(as.integer(fT$env$nAnalyticGradDirect), 0L)   # out of scope: no analytic gradient
  })

  test_that("generalized (Poisson) ll() fast fit matches the finite-difference fit", {
    skip_on_cran()
    set.seed(42)
    N <- 25L; nobs <- 6L
    sim <- do.call(rbind, lapply(seq_len(N), function(i) {
      x <- rnorm(nobs); e1 <- rnorm(1, 0, sqrt(0.4)); e2 <- rnorm(1, 0, sqrt(0.2))
      data.frame(ID = i, TIME = seq_len(nobs),
                 DV = rpois(nobs, exp(1.2 + e1 + (0.5 + e2) * x)), x = x, EVID = 0)
    }))
    pois <- function() {
      ini({ tint <- 1.2; tslp <- 0.5; eta.int ~ 0.4; eta.slp ~ 0.2 })
      model({ lam <- exp(tint + eta.int + (tslp + eta.slp) * x)
              ll(cp) ~ DV * log(lam) - lam - lgamma(DV + 1) })
    }
    # sigdig pinned at 4: at the default 3 this model's own optimizer noise is ~9% in
    # theta (measured), which is larger than the fast-vs-fd difference under test.
    .ctl <- function(fast) foceiControl(print = 0L, covMethod = "", sigdig = 4, fast = fast)
    ## cached fast=FALSE reference -- see helper-gradref.R
    fF <- .numRef("fit-fd-ll-pois", function() {
      .f <- suppressMessages(nlmixr2(pois, sim, "focei", .ctl(FALSE)))
      list(objf = .f$objf, theta = unname(.f$theta))
    })
    fT <- suppressMessages(nlmixr2(pois, sim, "focei", .ctl(TRUE)))
    expect_equal(as.numeric(fT$objf), as.numeric(fF$objf), tolerance = 1e-2)
    expect_equal(unname(fT$theta), fF$theta, tolerance = 2e-2)
    # ODE-free models pool too: the augmented model has no ODE states, but it has real
    # lhs outputs and its sensitivities are plain symbolic derivatives.
    expect_gt(fT$env$nAnalyticGradDirect, 0L)
  })

  ## ---- multiple endpoints (#838) -------------------------------------------------
  ## An ll() endpoint written as the exact normal log-density is the SAME likelihood as
  ## the equivalent add() endpoint, so the two objectives may differ only by the small
  ## log|H| term (FOCEi's Gauss-Newton Hessian vs the ll() path's exact one).  That makes
  ## the Gaussian twin a self-validating reference -- no checked-in constant to drift.
  ##
  ## It is what caught #838: the endpoint's distribution was read one row before
  ## calc_lhs() evaluated rx_yj_, so a subject's FIRST observation was scored as normal.
  ## For a general-likelihood endpoint that treats rx_pred_ (a log density, down to -557
  ## here) as a prediction against DV with rx_r_ = 0 forced to variance 1, contributing
  ## -0.5*(logDensity - DV)^2.  Measured objf 11,463,665.75 against the twin's 53,695.56.
  ## Single-endpoint models were immune: with one endpoint rx_yj_ is a constant that is
  ## already installed at init, so nothing depended on the stale read.
  .pkpdLL <- function() {
    ini({ tka <- 0.5; tcl <- -2; tv <- 2; emax <- 2; ec50 <- 1
          add.pk <- 1; add.pd <- 3; eta.cl ~ 0.1 })
    model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
            d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
            cp <- center / v; pca <- emax * cp / (ec50 + cp)
            ll(cp)  ~ -0.5 * log(2 * pi) - log(add.pk) - 0.5 * ((DV - cp) / add.pk)^2
            ll(pca) ~ -0.5 * log(2 * pi) - log(add.pd) - 0.5 * ((DV - pca) / add.pd)^2 })
  }
  .pkpdGauss <- function() {
    ini({ tka <- 0.5; tcl <- -2; tv <- 2; emax <- 2; ec50 <- 1
          add.pk <- 1; add.pd <- 3; eta.cl ~ 0.1 })
    model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
            d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
            cp <- center / v; pca <- emax * cp / (ec50 + cp)
            cp ~ add(add.pk) | cp
            pca ~ add(add.pd) | pca })
  }
  .phCtl <- function(fast) foceiControl(print = 0L, covMethod = "", fast = fast, sigdig = 4,
                                        maxOuterIterations = 0L, maxInnerIterations = 100L)

  test_that("multi-endpoint ll() objective matches its Gaussian twin (#838)", {
    skip_on_cran(); skip_if_not_installed("nlmixr2data")
    d <- nlmixr2data::warfarin
    ## The twin is fitted LIVE, deliberately not cached through helper-gradref.R.  It is
    ## the reference's whole value that it is derived, not frozen: a cached twin recorded
    ## during a bad run silently becomes the thing under test (observed once -- objf
    ## 61316 with every eta exactly 0 -- and it persists).
    ll <- suppressMessages(suppressWarnings(nlmixr2(.pkpdLL, d, "focei", .phCtl(TRUE))))
    gs <- suppressMessages(suppressWarnings(nlmixr2(.pkpdGauss, d, "focei", .phCtl(TRUE))))
    # guard the reference itself: a collapsed fit must fail loudly here, not silently
    # agree with a collapsed ll() arm
    expect_true(all(abs(gs$eta$eta.cl) > 1e-8))
    # the whole failure was a data-sized excess (1.1e7 on an objective of 5.4e4), so a
    # loose relative tolerance still pins it; log|H| differs by O(1) between the two
    expect_equal(as.numeric(ll$objf), as.numeric(gs$objf), tolerance = 1e-3)
    # per subject too -- a compensating error across subjects would pass the total
    expect_equal(ll$env$etaObf$OBJI, gs$env$etaObf$OBJI, tolerance = 1e-2)
    # and the conditional estimates agree (they differed by 0.17 while the bug was live)
    expect_equal(ll$eta$eta.cl, gs$eta$eta.cl, tolerance = 1e-2)
    # the analytic gradient really was used for the ll() arm
    expect_gt(ll$env$nAnalyticGradDirect, 0L)
  })

  test_that("multi-endpoint ll() analytic gradient matches FD (#838)", {
    skip_on_cran(); skip_if_not_installed("nlmixr2data")
    d <- nlmixr2data::warfarin
    ph <- suppressMessages(suppressWarnings(nlmixr2(.pkpdLL, d, "focei", .phCtl(TRUE))))
    g <- .foceiGradDirect(ph)
    expect_false(is.null(g))
    base <- fixef(ph)
    ofvAt <- function(nm, val) {
      ui2 <- do.call(rxode2::ini, c(list(ph$finalUi), setNames(list(val), nm)))
      suppressMessages(suppressWarnings(nlmixr2(ui2, d, "focei", .phCtl(TRUE))))$objf
    }
    h <- 1e-3
    ## cached: a property of the model/data/theta, not of the gradient -- helper-gradref.R
    fd <- .gradRef("ll-multiple-endpoint", function()
      vapply(names(base), function(nm) (ofvAt(nm, base[nm] + h) - ofvAt(nm, base[nm] - h)) / (2 * h),
             numeric(1)))
    # 2% -- the reference's own step noise is ~1% on tcl (h=1e-3 vs 1e-4); before the
    # objective fix the worst component was off by 373x
    expect_equal(unname(g[names(base)]), unname(fd), tolerance = 0.02)
  })

})
