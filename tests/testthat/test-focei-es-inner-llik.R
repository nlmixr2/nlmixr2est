nmTest({
  # An ll() fit registers innerLlik as the inner model, and focei.R's fit-wide
  # event-sensitivity install must bind that model, re-resolved on the library the
  # fit loaded: otherwise a compacted inner solve runs another model's dydt (a heap
  # overflow, "free(): invalid next size").  The dose's bioavailability depends on
  # an eta, so the jump sensitivities are live through the fit and the covariance.
  test_that("an ll() fit with eta-dependent bioavailability fits and has a covariance", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    mod <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45; tf <- 0
        eta.cl ~ 0.3; eta.f ~ 0.1
        lsd <- log(0.7)
      })
      model({
        ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        f(depot) <- exp(tf + eta.f)
        cp <- center / v
        ll(err) ~ -0.5 * log(2 * pi) - lsd - 0.5 * ((DV - cp) / exp(lsd))^2
      })
    }
    f <- suppressMessages(suppressWarnings(
      nlmixr2(mod, nlmixr2data::theo_sd, "focei", foceiControl(print = 0L, maxOuterIterations = 5L))
    ))
    expect_true(is.finite(f$objf))
    expect_true(all(is.finite(sqrt(diag(f$cov)))))
  })
})
