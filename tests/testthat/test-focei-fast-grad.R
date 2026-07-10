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
  # default outer optimizer is nlminb for both base and fast (fast=TRUE with
  # lbfgsb3c can settle at a worse optimum on ill-conditioned models)
  expect_equal(foceiControl()$outerOpt, -1L)                 # nlminb -> custom (-1)
  expect_equal(foceiControl(fast = TRUE)$outerOpt, -1L)      # nlminb -> custom (-1)
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

test_that("FOCE (nonmem) analytic gradient matches central differences", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("nlmixr2data")
  # add+prop so the frozen-R0 depends on the population prediction (exercises the
  # nonmem a0-chain) while staying > 0 (no vanishing-variance guard)
  offAP <- function() {
    ini({ tka <- 0.2; tcl <- 1.2; tv <- 3.2; eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
          add.sd <- 0.4; prop.sd <- 0.15 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
            d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
            cp <- center / v; cp ~ add(add.sd) + prop(prop.sd) })
  }
  d <- nlmixr2data::theo_sd
  ph <- suppressMessages(nlmixr2(offAP, d, "foce",
        foceiControl(print = 0L, covMethod = "", fast = TRUE,
                     maxOuterIterations = 0L, maxInnerIterations = 300L)))
  g <- .foceiGradAnalyticCalc(ph)
  expect_false(is.null(g))
  base <- fixef(ph)
  ofvAt <- function(nm, val) {
    ui2 <- do.call(rxode2::ini, c(list(ph$finalUi), setNames(list(val), nm)))
    suppressMessages(suppressWarnings(nlmixr2(ui2, d, "foce",
      foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L,
                   maxInnerIterations = 300L))))$objf
  }
  h <- 1e-3
  fd <- vapply(names(base), function(nm) (ofvAt(nm, base[nm] + h) - ofvAt(nm, base[nm] - h)) / (2 * h), numeric(1))
  expect_equal(unname(g[names(base)]), unname(fd), tolerance = 0.01)
})

test_that("FOCE censored (M3) analytic gradient matches central differences", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("nlmixr2data")
  # censored FOCE re-solves the EBE with the exact censored rho_f/rho_ff at the frozen R0
  # (.foceiAnalyticFoceEbe) and uses the censored score + R0-chain cross deriv rfR
  offAP <- function() {
    ini({ tka <- 0.2; tcl <- 1.2; tv <- 3.2; eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
          add.sd <- 0.4; prop.sd <- 0.15 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
            d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
            cp <- center / v; cp ~ add(add.sd) + prop(prop.sd) })
  }
  d <- nlmixr2data::theo_sd
  d$CENS <- ifelse(d$DV < 2 & d$EVID == 0, 1L, 0L); d$DV[d$CENS == 1] <- 2
  ph <- suppressMessages(nlmixr2(offAP, d, "foce",
        foceiControl(print = 0L, covMethod = "", fast = TRUE,
                     maxOuterIterations = 0L, maxInnerIterations = 300L)))
  g <- .foceiGradAnalyticCalc(ph)
  expect_false(is.null(g))
  base <- fixef(ph)
  ofvAt <- function(nm, val) {
    ui2 <- do.call(rxode2::ini, c(list(ph$finalUi), setNames(list(val), nm)))
    suppressMessages(suppressWarnings(nlmixr2(ui2, d, "foce",
      foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L,
                   maxInnerIterations = 300L))))$objf
  }
  h <- 1e-3
  fd <- vapply(names(base), function(nm) (ofvAt(nm, base[nm] + h) - ofvAt(nm, base[nm] - h)) / (2 * h), numeric(1))
  expect_equal(unname(g[names(base)]), unname(fd), tolerance = 0.02)
})

test_that("estimated boxCox/yeoJohnson lambda: analytic gradient matches central differences", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("nlmixr2data")
  # both-sides transform with an ESTIMATED lambda: lambda is a theta-like direction
  # (df'/dlambda) plus the DV-transform residual chain (dy'/dlambda) and the -2 log|J|
  # Jacobian.  Off initials so every gradient (incl. lambda) is large-signal.
  d <- nlmixr2data::theo_sd
  mkBox <- function() {
    ini({ tka <- 0.15; tcl <- 1.25; tv <- 3.15; eta.ka ~ 0.5; eta.cl ~ 0.25; eta.v ~ 0.1
          add.sd <- 0.75; lambda <- c(-1, 0.7, 2) })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
            d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
            cp <- center / v; cp ~ add(add.sd) + boxCox(lambda) })
  }
  mkYj <- function() {
    ini({ tka <- 0.15; tcl <- 1.25; tv <- 3.15; eta.ka ~ 0.5; eta.cl ~ 0.25; eta.v ~ 0.1
          add.sd <- 0.75; lambda <- c(-1, 0.7, 2) })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
            d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
            cp <- center / v; cp ~ add(add.sd) + yeoJohnson(lambda) })
  }
  chk <- function(mk, est) {
    ph <- suppressMessages(nlmixr2(mk, d, est,
          foceiControl(print = 0L, covMethod = "", fast = TRUE,
                       maxOuterIterations = 0L, maxInnerIterations = 300L)))
    g <- .foceiGradAnalyticCalc(ph)
    expect_false(is.null(g))
    base <- fixef(ph)
    ofvAt <- function(nm, val) {
      ui2 <- do.call(rxode2::ini, c(list(ph$finalUi), setNames(list(val), nm)))
      suppressMessages(suppressWarnings(nlmixr2(ui2, d, est,
        foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L,
                     maxInnerIterations = 300L))))$objf
    }
    h <- 1e-3
    fd <- vapply(names(base), function(nm) (ofvAt(nm, base[nm] + h) - ofvAt(nm, base[nm] - h)) / (2 * h), numeric(1))
    expect_equal(unname(g[names(base)]), unname(fd), tolerance = 0.01)
  }
  chk(mkBox, "focei"); chk(mkYj, "focei"); chk(mkYj, "foce")
})

test_that("analytic outer gradient matches FD for a covariate model", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("nlmixr2data")
  # a covariate (wtCl*WT) in the structural model: exercises the covariate direction
  # and the param() covariate declaration in the augmented outer model
  set.seed(1)
  d <- do.call(rbind, lapply(1:12, function(i)
    data.frame(ID = i, TIME = c(0, .5, 1, 2, 4, 8), EVID = c(101, 0, 0, 0, 0, 0),
               AMT = c(100, 0, 0, 0, 0, 0), DV = c(NA, 8, 9, 7, 4, 1) + rnorm(6, 0, .3),
               WT = runif(1, 50, 90))))
  covm <- function() {
    ini({ tka <- 0.2; tcl <- 1.2; tv <- 3.2; wtCl <- 0.01; eta.cl ~ 0.3; prop.sd <- 0.2 })
    model({ ka <- exp(tka); cl <- exp(tcl + eta.cl + wtCl * WT); v <- exp(tv)
            d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
            cp <- center / v; cp ~ prop(prop.sd) })
  }
  ph <- suppressMessages(suppressWarnings(nlmixr2(covm, d, "focei",
        foceiControl(print = 0L, covMethod = "", fast = TRUE,
                     maxOuterIterations = 0L, maxInnerIterations = 200L))))
  g <- .foceiGradAnalyticCalc(ph)
  expect_false(is.null(g))
  base <- fixef(ph)
  ofvAt <- function(nm, val) {
    ui2 <- do.call(rxode2::ini, c(list(ph$finalUi), setNames(list(val), nm)))
    suppressMessages(suppressWarnings(nlmixr2(ui2, d, "focei",
      foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L,
                   maxInnerIterations = 200L))))$objf
  }
  h <- 1e-3
  fd <- vapply(names(base), function(nm) (ofvAt(nm, base[nm] + h) - ofvAt(nm, base[nm] - h)) / (2 * h), numeric(1))
  expect_equal(unname(g[names(base)]), unname(fd), tolerance = 0.01)
})

test_that("analytic outer gradient matches FD for a multiple-endpoint model", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("nlmixr2data")
  # two modeled endpoints (PK cp + PD pca): rx_pred_/rx_r_ are single dvid-conditional
  # expressions, so solving against the dataset selects each endpoint's prediction and
  # variance per observation -- the (f,R) path handles both endpoints' sigmas
  d <- nlmixr2data::warfarin
  pkpd <- function() {
    ini({ tka <- 0.5; tcl <- -2; tv <- 2; emax <- 2; ec50 <- 1; add.pk <- 1; add.pd <- 3; eta.cl ~ 0.1 })
    model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
            d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
            cp <- center / v; pca <- emax * cp / (ec50 + cp)
            cp ~ add(add.pk) | cp
            pca ~ add(add.pd) | pca })
  }
  ph <- suppressMessages(suppressWarnings(nlmixr2(pkpd, d, "focei",
        foceiControl(print = 0L, covMethod = "", fast = TRUE,
                     maxOuterIterations = 0L, maxInnerIterations = 100L))))
  g <- .foceiGradAnalyticCalc(ph)
  expect_false(is.null(g))
  base <- fixef(ph)
  ofvAt <- function(nm, val) {
    ui2 <- do.call(rxode2::ini, c(list(ph$finalUi), setNames(list(val), nm)))
    suppressMessages(suppressWarnings(nlmixr2(ui2, d, "focei",
      foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L,
                   maxInnerIterations = 100L))))$objf
  }
  h <- 1e-3
  fd <- vapply(names(base), function(nm) (ofvAt(nm, base[nm] + h) - ofvAt(nm, base[nm] - h)) / (2 * h), numeric(1))
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

test_that("modeled dosing parameters (f/lag) use jump sensitivities and match FD", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("nlmixr2data")
  d <- nlmixr2data::theo_sd
  # bioavailability f() and absorption lag() modeled on theta/eta: the outer gradient
  # (and its EBE derivative) needs the dose-based second-order jump sensitivities.
  mDose <- function() {
    ini({ tka <- 0.3; tcl <- 1.1; tv <- 3.3; tf <- 0.1; tl <- -1.5
          eta.ka ~ 0.4; eta.cl ~ 0.2; eta.f ~ 0.1; add.sd <- 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
            fdepot <- exp(tf + eta.f); ld <- exp(tl)
            d/dt(depot) <- -ka * depot; f(depot) <- fdepot; lag(depot) <- ld
            d/dt(center) <- ka * depot - cl / v * center; cp <- center / v; cp ~ add(add.sd) })
  }
  ph <- suppressMessages(suppressWarnings(nlmixr2(mDose, d, "focei",
        foceiControl(print = 0L, covMethod = "", fast = TRUE,
                     maxOuterIterations = 0L, maxInnerIterations = 300L))))
  g <- .foceiGradAnalyticCalc(ph)
  expect_false(is.null(g))                                    # jump sensitivities keep it in scope
  base <- fixef(ph)
  ofvAt <- function(nm, val) {
    ui2 <- do.call(rxode2::ini, c(list(ph$finalUi), setNames(list(val), nm)))
    suppressMessages(suppressWarnings(nlmixr2(ui2, d, "focei",
      foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L,
                   maxInnerIterations = 300L))))$objf
  }
  h <- 1e-3
  fd <- vapply(names(base), function(nm) (ofvAt(nm, base[nm] + h) - ofvAt(nm, base[nm] - h)) / (2 * h), numeric(1))
  expect_equal(unname(g[names(base)]), unname(fd), tolerance = 0.02)
})

test_that("mceta=-2 (Eq-48) is the default and all fast mceta modes agree", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("nlmixr2data")
  expect_equal(foceiControl()$mceta, -2L)                    # new global default
  d <- nlmixr2data::theo_sd
  ofv <- vapply(c(-2L, -1L, 0L), function(mc)
    suppressMessages(nlmixr2(.fast_one_cmt, d, "focei",
      foceiControl(print = 0L, covMethod = "", fast = TRUE, mceta = mc)))$objf, numeric(1))
  # Eq-48 extrapolation / jump / reset must all reach the same optimum
  expect_equal(ofv[1], ofv[2], tolerance = 0.02)
  expect_equal(ofv[1], ofv[3], tolerance = 0.02)
  # mceta=-2 with fast=FALSE degrades to keep-last-eta (no analytic sensitivity) and fits
  f <- suppressMessages(nlmixr2(.fast_one_cmt, d, "focei",
        foceiControl(print = 0L, covMethod = "", fast = FALSE, mceta = -2L)))
  expect_true(is.finite(f$objf))
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
