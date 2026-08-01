# Analytic FOCEI outer gradient (foceiControl(fast=TRUE)): the analytic gradient
# agrees with central differences of the objective, a fast fit matches a
# finite-difference fit, out-of-scope models fall back transparently, and the
# fast/derivative-free control defaults behave.
#
# Weekly-batched via .slowBatches in tests/testthat.R -- do NOT add skip_on_ci().
# The weekly runner also sets CI=true, so skip_on_ci() here would skip these
# everywhere and leave the fast/analytic-gradient path with no CI coverage.

nmTest({
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
    # default outer optimizer: nlminb for finite differences, lbfgsb3c for the
    # analytic ("fast") gradient
    expect_equal(foceiControl()$outerOpt, -1L)                 # nlminb -> custom (-1)
    expect_equal(foceiControl(fast = TRUE)$outerOpt, 1L)       # lbfgsb3c
    expect_equal(foceiControl(fast = TRUE)$outerOptTxt, "lbfgsb3c")
    expect_true(foceiControl(fast = TRUE)$fast)
    expect_false(foceiControl()$fast)
    # an explicit outerOpt still wins under fast
    expect_equal(foceiControl(fast = TRUE, outerOpt = "nlminb")$outerOpt, -1L)
    # a defaulted optimizer re-defaults under a *f wrapper; an explicit one is kept
    expect_equal(nlmixr2est:::.foceiFastCtl(list(foceiControl()), foceiControl)$outerOptTxt,
                 "lbfgsb3c")
    expect_equal(nlmixr2est:::.foceiFastCtl(list(foceiControl(outerOpt = "nlminb")), foceiControl)$outerOptTxt,
                 "nlminb")
    # derivative-free outerOpt + fast -> fast cleared with a warning
    expect_warning(.c <- foceiControl(fast = TRUE, outerOpt = "bobyqa"), "derivative-free")
    expect_false(.c$fast)
  })

  test_that("analytic outer gradient matches central differences (theta + sigma)", {
    skip_on_cran()
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
          foceiControl(print = 0L, covMethod = "", fast = TRUE, sigdig = 4,
                       maxOuterIterations = 0L, maxInnerIterations = 300L)))
    g <- .foceiGradDirect(ph)
    expect_false(is.null(g))
    base <- fixef(ph)
    ofvAt <- function(nm, val) {
      ui2 <- do.call(rxode2::ini, c(list(ph$finalUi), setNames(list(val), nm)))
      suppressMessages(suppressWarnings(nlmixr2(ui2, d, "focei",
        foceiControl(print = 0L, covMethod = "", sigdig = 4, maxOuterIterations = 0L,
                     maxInnerIterations = 300L))))$objf
    }
    h <- 1e-3
    ## cached reference -- see helper-gradref.R
    fd <- .gradRef("subfinal-theta-sigma", function()
      vapply(names(base), function(nm) (ofvAt(nm, base[nm] + h) - ofvAt(nm, base[nm] - h)) / (2 * h), numeric(1)))
    # large-signal gradients: analytic vs central-difference within 1% relative
    expect_equal(unname(g[names(base)]), unname(fd), tolerance = 0.01)
  })

  test_that("estimated boxCox/yeoJohnson lambda: analytic gradient matches central differences", {
    skip_on_cran()
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
    chk <- function(mk, est, nm) {
      ph <- suppressMessages(nlmixr2(mk, d, est,
            foceiControl(print = 0L, covMethod = "", fast = TRUE, sigdig = 4,
                         maxOuterIterations = 0L, maxInnerIterations = 300L)))
      g <- .foceiGradDirect(ph)
      expect_false(is.null(g))
      base <- fixef(ph)
      ofvAt <- function(nm, val) {
        ui2 <- do.call(rxode2::ini, c(list(ph$finalUi), setNames(list(val), nm)))
        suppressMessages(suppressWarnings(nlmixr2(ui2, d, est,
          foceiControl(print = 0L, covMethod = "", sigdig = 4, maxOuterIterations = 0L,
                       maxInnerIterations = 300L))))$objf
      }
      h <- 1e-3
      ## cached reference -- see helper-gradref.R
      fd <- .gradRef(paste0("subfinal-lambda-", nm), function()
        vapply(names(base), function(nm) (ofvAt(nm, base[nm] + h) - ofvAt(nm, base[nm] - h)) / (2 * h), numeric(1)))
      expect_equal(unname(g[names(base)]), unname(fd), tolerance = 0.01)
    }
    ## chk(mkYj, "foce") removed: FOCE is out of scope for the analytic gradient (#836)
    chk(mkBox, "focei", "boxcox"); chk(mkYj, "focei", "yeojohnson")
  })

  test_that("analytic outer gradient matches FD for a covariate model", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    # a covariate (wtCl*WT) in the structural model: exercises the covariate direction
    # and the param() covariate declaration in the augmented outer model
    .testSeed(1)
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
          foceiControl(print = 0L, covMethod = "", fast = TRUE, sigdig = 4,
                       maxOuterIterations = 0L, maxInnerIterations = 200L))))
    g <- .foceiGradDirect(ph)
    expect_false(is.null(g))
    base <- fixef(ph)
    ofvAt <- function(nm, val) {
      ui2 <- do.call(rxode2::ini, c(list(ph$finalUi), setNames(list(val), nm)))
      suppressMessages(suppressWarnings(nlmixr2(ui2, d, "focei",
        foceiControl(print = 0L, covMethod = "", sigdig = 4, maxOuterIterations = 0L,
                     maxInnerIterations = 200L))))$objf
    }
    h <- 1e-3
    ## cached reference -- see helper-gradref.R
    fd <- .gradRef("subfinal-covariate", function()
      vapply(names(base), function(nm) (ofvAt(nm, base[nm] + h) - ofvAt(nm, base[nm] - h)) / (2 * h), numeric(1)))
    expect_equal(unname(g[names(base)]), unname(fd), tolerance = 0.01)
  })

  test_that("analytic outer gradient matches FD for a multiple-endpoint model", {
    skip_on_cran()
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
          foceiControl(print = 0L, covMethod = "", fast = TRUE, sigdig = 4,
                       maxOuterIterations = 0L, maxInnerIterations = 100L))))
    g <- .foceiGradDirect(ph)
    expect_false(is.null(g))
    base <- fixef(ph)
    ofvAt <- function(nm, val) {
      ui2 <- do.call(rxode2::ini, c(list(ph$finalUi), setNames(list(val), nm)))
      suppressMessages(suppressWarnings(nlmixr2(ui2, d, "focei",
        foceiControl(print = 0L, covMethod = "", sigdig = 4, maxOuterIterations = 0L,
                     maxInnerIterations = 100L))))$objf
    }
    h <- 1e-3
    ## cached reference -- see helper-gradref.R
    fd <- .gradRef("subfinal-multiple-endpoint", function()
      vapply(names(base), function(nm) (ofvAt(nm, base[nm] + h) - ofvAt(nm, base[nm] - h)) / (2 * h), numeric(1)))
    expect_equal(unname(g[names(base)]), unname(fd), tolerance = 0.01)
  })

  test_that("fast=TRUE fit matches the finite-difference fit", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    d <- nlmixr2data::theo_sd
    ## cached fast=FALSE reference -- see helper-gradref.R
    f0 <- .numRef("fit-fd-one-cmt-focei-subfinal", function() {
      .f <- suppressMessages(nlmixr2(.fast_one_cmt, d, "focei",
              foceiControl(print = 0L, covMethod = "", fast = FALSE)))
      list(objf = .f$objf, fixef = unname(fixef(.f)))
    })
    fF <- suppressMessages(nlmixr2(.fast_one_cmt, d, "focei", foceiControl(print = 0L, covMethod = "", fast = TRUE)))
    expect_equal(fF$objf, f0$objf, tolerance = 0.02)
    expect_equal(unname(fixef(fF)), f0$fixef, tolerance = 1e-2)
    # the analytic gradient must actually be CONSUMED by the optimizer (a fit that
    # silently falls back to FD also "matches", so assert usage directly)
    .gt <- fF$parHistData$type
    expect_gt(sum(.gt == "Analytic Gradient"), 0)
    expect_equal(sum(.gt %in% c("Gill83 Gradient", "Mixed Gradient",
                                "Forward Difference", "Central Difference")), 0)
    expect_match(fF$extra, "grad: analytic")
  })

})
