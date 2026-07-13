## Phase 3: VAE decoder -- structural prediction f, residual variance R, and
## eta-sensitivities df/deta, dR/deta from the reused focei analytic sensitivity
## model, plus the ELBO p(x|z) data term and its eta-gradient. Validated on the
## theophylline closed-form 1-cmt solution and finite differences.

nmTest({
  test_that("vae decoder f/R/df-deta/dR-deta match closed form + FD (theophylline)", {
    theo <- function() {
      ini({
        lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03
        add.err <- 0.7
      })
      model({
        ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot
        d/dt(central) = ka * depot - ke * central
        cp <- central / V
        cp ~ add(add.err)
      })
    }
    ui <- rxode2::assertRxUi(theo)
    am <- .vaeDecoderModel(ui)
    expect_false(is.null(am))

    lka <- log(1.8); lke <- log(0.086); lV <- log(32); addErr <- 0.7
    th <- c(THETA_1_ = lka, THETA_2_ = lke, THETA_3_ = lV, THETA_4_ = addErr)
    eta <- c(0.15, -0.08, 0.2)
    dose <- 320
    times <- seq(0.25, 24, length.out = 12)
    ev <- rxode2::et(amt = dose, cmt = "depot", time = 0)
    ev <- rxode2::et(ev, times)

    E <- .vaeDecoderSolveSubject(am, th, eta, ev, times)
    expect_false(is.null(E))
    expect_equal(length(E$f), length(times))

    ## closed-form cp and its eta-derivatives
    cf <- function(eta, t) {
      ka <- exp(lka + eta[1]); ke <- exp(lke + eta[2]); V <- exp(lV + eta[3])
      dose * ka / (V * (ka - ke)) * (exp(-ke * t) - exp(-ka * t))
    }
    fCf <- cf(eta, times)
    expect_lt(max(abs(E$f - fCf)) / max(abs(fCf)), 1e-6)

    ## df/deta vs FD of the closed form
    h <- 1e-6
    for (k in 1:3) {
      ep <- eta; ep[k] <- ep[k] + h; em <- eta; em[k] <- em[k] - h
      fd <- (cf(ep, times) - cf(em, times)) / (2 * h)
      expect_lt(max(abs(E$a[, k] - fd)) / max(abs(fd)), 1e-4)
    }

    ## additive error: R = add.err^2, dR/deta = 0
    expect_lt(max(abs(E$R - addErr^2)), 1e-8)
    expect_lt(max(abs(E$aR)), 1e-8)

    ## ELBO data term p(x|z) gradient vs FD (synthetic DV = closed form + noise)
    set.seed(42)
    y <- fCf + rnorm(length(fCf), 0, 0.3)
    px <- .vaeDecoderPxz(E, y)
    pxzOf <- function(eta) {
      Ee <- .vaeDecoderSolveSubject(am, th, eta, ev, times)
      .vaeDecoderPxz(Ee, y)$pxz
    }
    for (k in 1:3) {
      ep <- eta; ep[k] <- ep[k] + h; em <- eta; em[k] <- em[k] - h
      fd <- (pxzOf(ep) - pxzOf(em)) / (2 * h)
      expect_lt(abs(px$gEta[k] - fd) / max(abs(fd), 1), 1e-4)
    }
  })
})
