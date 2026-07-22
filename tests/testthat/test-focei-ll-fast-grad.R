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
    is <- function(m) nlmixr2est:::.foceiLLGradInScope(rxode2::rxUiDecompress(m()))
    expect_true(is(.ll_ode))                 # ODE log-likelihood endpoint -> in scope
    expect_false(is(.gauss_ode))             # Gaussian -> the (f,R) analytic path, not the ll path
    expect_true(is(.ll_lincmt))              # linCmt() passes the coarse gate (falls back to FD at build)
  })

})
