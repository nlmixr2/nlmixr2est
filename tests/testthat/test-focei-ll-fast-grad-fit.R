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
    fF <- suppressMessages(nlmixr2(.ll_ode, d, "focei", foceiControl(print = 0L, covMethod = "", fast = FALSE)))
    fT <- suppressMessages(nlmixr2(.ll_ode, d, "focei", foceiControl(print = 0L, covMethod = "", fast = TRUE)))
    # fast=TRUE for an ll() endpoint uses the exact-Hessian objective + analytic outer
    # gradient; it converges to the same MLE within optimizer tolerance
    expect_equal(as.numeric(fT$objf), as.numeric(fF$objf), tolerance = 1e-2)
    expect_equal(unname(fT$theta), unname(fF$theta), tolerance = 1e-2)
  })

  test_that("ll() with linCmt() falls back gracefully to finite differences", {
    skip_on_cran(); skip_if_not_installed("nlmixr2data")
    # linCmt() has no 2nd-order state sensitivities (solved form), so both the exact-Hessian
    # objective and the analytic outer gradient are out of scope; fast=TRUE must transparently
    # fall back to finite differences and still converge to the same MLE (not error).
    d <- nlmixr2data::theo_sd
    fF <- suppressMessages(nlmixr2(.ll_lincmt, d, "focei", foceiControl(print = 0L, covMethod = "", fast = FALSE)))
    fT <- suppressMessages(nlmixr2(.ll_lincmt, d, "focei", foceiControl(print = 0L, covMethod = "", fast = TRUE)))
    expect_equal(as.numeric(fT$objf), as.numeric(fF$objf), tolerance = 1e-2)
    expect_equal(unname(fT$theta), unname(fF$theta), tolerance = 1e-2)
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
    fF <- suppressMessages(nlmixr2(pois, sim, "focei", foceiControl(print = 0L, covMethod = "", fast = FALSE)))
    fT <- suppressMessages(nlmixr2(pois, sim, "focei", foceiControl(print = 0L, covMethod = "", fast = TRUE)))
    expect_equal(as.numeric(fT$objf), as.numeric(fF$objf), tolerance = 1e-2)
    expect_equal(unname(fT$theta), unname(fF$theta), tolerance = 2e-2)
  })

})
