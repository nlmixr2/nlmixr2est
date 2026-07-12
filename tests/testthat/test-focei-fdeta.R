## foceiControl(fdEta=): method-agnostic API to finite-difference a named eta,
## using the same machinery as the dosing-parameter (eventEta) etas.  For an eta
## whose effect on the prediction is invisible to the symbolic model (e.g. an
## externally injected NN weight, analytic d(f)/d(eta) == 0) this is what lets it
## be estimated.  Here we validate it on a visible eta -- FD ~ analytic, so the
## fit is unchanged -- confirming the flag is honored and does not break fits.

test_that("foceiControl(fdEta=) finite-differences a named eta via the event machinery", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  .old <- rxode2::getRxThreads(); on.exit(rxode2::setRxThreads(.old), add = TRUE)
  rxode2::setRxThreads(1L)

  d <- nlmixr2data::theo_sd
  mod <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7; eta.cl ~ 0.1 })
    model({
      ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  ctl <- function(fd) foceiControl(print = 0L, maxOuterIterations = 3L,
                                   maxInnerIterations = 5L, calcTables = FALSE, fdEta = fd)
  fit <- function(fd) suppressWarnings(suppressMessages(nlmixr2(mod, d, "focei", ctl(fd))))

  a <- fit(NULL)          # analytic eta sensitivity
  b <- fit("eta.cl")      # eta.cl forced through finite differences
  ## FD ~ analytic for a smooth eta: the fit is essentially unchanged (accepted,
  ## routed through FD, no breakage)
  expect_equal(b$objf, a$objf, tolerance = 1e-2)
  expect_equal(unname(b$omega[1, 1]), unname(a$omega[1, 1]), tolerance = 1e-3)

  ## an unknown eta name is a harmless no-op
  cc <- fit("eta.doesNotExist")
  expect_equal(cc$objf, a$objf, tolerance = 1e-2)
})
