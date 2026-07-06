# Analytic FOCEI outer gradient (foceiControl(fast=TRUE)): the analytic gradient
# agrees with central differences of the objective, a fast fit matches a
# finite-difference fit, out-of-scope models fall back transparently, and the
# fast/derivative-free control defaults behave.

.fast_one_cmt <- function() {
  ini({
    tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
    eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

test_that("fast=TRUE control defaults: outerOpt + derivative-free downgrade", {
  # base default outer optimizer is nlminb; fast default is lbfgsb3c
  expect_equal(foceiControl()$outerOpt, -1L)                 # nlminb -> custom (-1)
  expect_equal(foceiControl(fast = TRUE)$outerOpt, 1L)       # lbfgsb3c
  expect_true(foceiControl(fast = TRUE)$fast)
  expect_false(foceiControl()$fast)
  # an explicit outerOpt still wins under fast
  expect_equal(foceiControl(fast = TRUE, outerOpt = "nlminb")$outerOpt, -1L)
  # derivative-free outerOpt + fast -> fast cleared with a warning
  expect_warning(.c <- foceiControl(fast = TRUE, outerOpt = "bobyqa"), "derivative-free")
  expect_false(.c$fast)
})

test_that("analytic outer gradient matches central differences (theta + sigma)", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("nlmixr2data")
  # posthoc (eta*-only) at deliberately off initials so gradients are large-signal
  off <- function() {
    ini({ tka <- 0.2; tcl <- 1.2; tv <- 3.2; eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.9 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
            d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
            cp <- center / v; cp ~ add(add.sd) })
  }
  d <- nlmixr2data::theo_sd
  ph <- suppressMessages(nlmixr2(off, d, "focei",
        foceiControl(print = 0L, covMethod = "", fast = TRUE,
                     maxOuterIterations = 0L, maxInnerIterations = 300L)))
  g <- .foceiGradAnalyticCalc(ph)
  expect_false(is.null(g))
  base <- fixef(ph)
  ofvAt <- function(nm, val) {
    ui2 <- do.call(rxode2::ini, c(list(ph$finalUi), setNames(list(val), nm)))
    suppressMessages(suppressWarnings(nlmixr2(ui2, d, "focei",
      foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L,
                   maxInnerIterations = 300L))))$objf
  }
  h <- 1e-3
  fd <- vapply(names(base), function(nm) (ofvAt(nm, base[nm] + h) - ofvAt(nm, base[nm] - h)) / (2 * h), numeric(1))
  # large-signal gradients: analytic vs central-difference within 1% relative
  expect_equal(unname(g[names(base)]), unname(fd), tolerance = 0.01)
})

test_that("fast=TRUE fit matches the finite-difference fit", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("nlmixr2data")
  d <- nlmixr2data::theo_sd
  f0 <- suppressMessages(nlmixr2(.fast_one_cmt, d, "focei", foceiControl(print = 0L, covMethod = "", fast = FALSE)))
  fF <- suppressMessages(nlmixr2(.fast_one_cmt, d, "focei", foceiControl(print = 0L, covMethod = "", fast = TRUE)))
  expect_equal(fF$objf, f0$objf, tolerance = 0.02)
  expect_equal(unname(fixef(fF)), unname(fixef(f0)), tolerance = 1e-2)
})

test_that("out-of-scope model (linCmt) falls back to the finite-difference gradient", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  lin <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); linCmt() ~ add(add.sd) })
  }
  d <- nlmixr2data::theo_sd
  fF <- suppressMessages(suppressWarnings(nlmixr2(lin, d, "focei",
        foceiControl(print = 0L, covMethod = "", fast = TRUE))))
  expect_true(is.finite(fF$objf))                            # completes via FD fallback
})
