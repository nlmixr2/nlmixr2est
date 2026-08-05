# Scope gate for the analytic FOCEi outer gradient (foceiControl(fast=TRUE)) on general
# log-likelihood (ll()) / generalized endpoints.  The end-to-end fits that exercise the
# analytic objective + gradient live in the weekly-batched test-focei-ll-fast-grad-fit.R.

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
  .gauss_ode <- function() {
    ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
          eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
            d/dt(depot)  <- -ka * depot
            d/dt(center) <-  ka * depot - cl / v * center
            cp <- center / v; cp ~ add(add.sd) })
  }
  .ll_lincmt <- function() {
    ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
          eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
            cp <- linCmt()
            ll(err) ~ -0.5 * log(2 * pi) - log(add.sd) - 0.5 * ((DV - cp) / add.sd)^2 })
  }

  test_that("ll() analytic-gradient scope gates", {
    is <- function(m) .foceiLLGradInScope(rxode2::rxUiDecompress(rxode2::rxode2(m)))
    expect_true(is(.ll_ode))                 # ODE log-likelihood endpoint -> in scope
    expect_false(is(.gauss_ode))             # Gaussian -> the (f,R) analytic path, not the ll path
    expect_true(is(.ll_lincmt))              # linCmt() passes the coarse gate (falls back to FD at build)
  })

  test_that("the augmented sensitivity model builds for an ODE-free model", {
    # A purely algebraic model (no d/dt()) has NO state sensitivities and needs none:
    # d(rx_pred_)/ddir is a plain symbolic derivative.  The builder used to bail on the
    # empty .rxSens expansion, which sent every such model to finite differences.
    .algebraic <- function() {
      ini({ tint <- 1.2; tslp <- 0.5; eta.int ~ 0.4; eta.slp ~ 0.2 })
      model({ lam <- exp(tint + eta.int + (tslp + eta.slp) * x)
              ll(cp) ~ DV * log(lam) - lam - lgamma(DV + 1) })
    }
    .ui <- suppressWarnings(rxode2::rxUiDecompress(rxode2::rxode2(.algebraic)))
    expect_equal(length(rxode2::rxStateOde(.ui$loadPruneSens)), 0L)   # really has no states
    .d <- .foceiOuterDirsLL(.ui)
    expect_false(is.null(.d))
    .am <- .foceiAnalyticAugModelDirs(.ui, .d$dirs)
    expect_false(is.null(.am))
    expect_true(inherits(.am$augMod, "rxode2"))
    expect_equal(.am$ndir, length(.d$dirs))
    # and the reader can map its columns -- the pooled C++ solve reads through this
    expect_false(is.null(.vaeOuterCols(.am)))
  })

  test_that("multiple endpoints are in ll() analytic-gradient scope (#838)", {
    # Gated off while the multi-endpoint ll() OBJECTIVE was wrong -- the endpoint's
    # distribution was read a row early, so one observation per subject was scored as
    # normal.  With that fixed the gradient verifies, so the gate is lifted; the fits
    # that check the numbers live in test-focei-ll-fast-grad-fit.R.
    .twoLL <- function() {
      ini({ tka <- 0.5; tcl <- -2; tv <- 2; emax <- 2; ec50 <- 1
            add.pk <- 1; add.pd <- 3; eta.cl ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
              d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
              cp <- center / v; pca <- emax * cp / (ec50 + cp)
              ll(cp)  ~ -0.5 * log(2 * pi) - log(add.pk) - 0.5 * ((DV - cp) / add.pk)^2
              ll(pca) ~ -0.5 * log(2 * pi) - log(add.pd) - 0.5 * ((DV - pca) / add.pd)^2 })
    }
    # a Gaussian PK endpoint alongside an ll() PD endpoint: still the ll() path, because
    # one non-normal endpoint makes the WHOLE objective a general likelihood
    .mixed <- function() {
      ini({ tka <- 0.5; tcl <- -2; tv <- 2; emax <- 2; ec50 <- 1
            add.pk <- 1; add.pd <- 3; eta.cl ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
              d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
              cp <- center / v; pca <- emax * cp / (ec50 + cp)
              cp ~ add(add.pk) | cp
              ll(pca) ~ -0.5 * log(2 * pi) - log(add.pd) - 0.5 * ((DV - pca) / add.pd)^2 })
    }
    .is <- function(m) .foceiLLGradInScope(rxode2::rxUiDecompress(rxode2::rxode2(m)))
    expect_true(suppressWarnings(.is(.twoLL)))
    expect_true(suppressWarnings(.is(.mixed)))
    # an all-Gaussian multi-endpoint model still belongs to the (f,R) path, not this one
    .twoGauss <- function() {
      ini({ tka <- 0.5; tcl <- -2; tv <- 2; emax <- 2; ec50 <- 1
            add.pk <- 1; add.pd <- 3; eta.cl ~ 0.1 })
      model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
              d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
              cp <- center / v; pca <- emax * cp / (ec50 + cp)
              cp ~ add(add.pk) | cp
              pca ~ add(add.pd) | pca })
    }
    expect_false(suppressWarnings(.is(.twoGauss)))
  })

})
