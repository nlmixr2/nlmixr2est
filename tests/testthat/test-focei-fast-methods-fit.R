# Fit-based checks for the "*f" convenience methods (each is its base method
# with fast=TRUE default): the analytic gradient is actually consumed and the
# fast fit matches the finite-difference fit.
#
# Split out of test-focei-fast-methods.R (which keeps the cheap registration /
# control checks in the push/PR subset).  These are full fits, so they are
# weekly-batched via .slowBatches in tests/testthat.R -- do NOT add
# skip_on_ci().  The weekly runner also sets CI=true, so skip_on_ci() here
# would skip them everywhere and leave the fast path with no CI coverage.

nmTest({
  .fastm_one_cmt <- function() {
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

  test_that("mu-referenced families use the analytic gradient (fast matches FD)", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    d <- nlmixr2data::theo_sd
    # mu-referenced covariate coefficient cl.wt (regression-updated); the analytic
    # gradient must exclude it from the outer optimizer's parameter set (A1b).
    mc <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; cl.wt <- 0.01; eta.ka ~ 0.3; eta.cl ~ 0.1; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl + cl.wt * (WT - 70)); v <- exp(tv)
              d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
              cp <- center / v; cp ~ add(add.sd) })
    }
    for (est in c("mfocei", "ifocei", "mfoce")) {
      ## cached fast=FALSE reference, per est -- see helper-gradref.R
      f0 <- .numRef(paste0("fit-fd-methods-", est), function()
        list(objf = suppressMessages(suppressWarnings(nlmixr2(mc, d, est,
               foceiControl(print = 0L, covMethod = "", fast = FALSE))))$objf))
      fF <- suppressMessages(suppressWarnings(nlmixr2(mc, d, est, foceiControl(print = 0L, covMethod = "", fast = TRUE))))
      expect_equal(fF$objf, f0$objf, tolerance = 0.05, info = est)
      # analytic gradient actually consumed on the mu-profiled parameter set
      .gt <- fF$parHistData$type
      expect_gt(sum(.gt == "Analytic Gradient"), 0)
      .nFd <- sum(.gt %in% c("Gill83 Gradient", "Mixed Gradient",
                             "Forward Difference", "Central Difference"))
      expect_match(fF$extra, "grad: analytic", info = est)
      if (est == "mfoce") {
        # FOCE may decline the analytic gradient on some evaluations -- measured 2 of
        # 7 here, BOTH attributed to the inner Newton solve (nGradDecline["newton"]).
        # What matters for issue #838 is not that a fallback never happens, but that
        # it is never SILENT: every decline must be attributed and the fit must say
        # so.  Assert that contract rather than a zero it does not promise.
        expect_gt(fF$env$nAnalyticGradDirect, 0, label = est)
        expect_equal(unname(fF$env$nFDGradFast), .nFd, info = est)
        expect_equal(sum(fF$env$nGradDecline), unname(fF$env$nFDGradFast), info = est)
        if (.nFd > 0) expect_match(fF$extra, "grad: analytic\\+fd", info = est)
        # See nlmixr2est issue on the FOCE Newton decline: if that is fixed this
        # method should reach a pure analytic gradient and join the branch below.
      } else {
        # mfocei / ifocei reach a PURE analytic gradient -- no fallback at all
        expect_equal(.nFd, 0, info = est)
        expect_equal(unname(fF$env$nFDGradFast), 0L, info = est)
        expect_false(grepl("analytic\\+fd", fF$extra), label = est)
      }
      expect_match(fF$extra, if (grepl("^i", est)) "mu: irls" else "mu: lin", info = est)
    }
  })

  test_that("est='foceif' equals est='focei' + fast=TRUE", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    d <- nlmixr2data::theo_sd
    fF   <- suppressMessages(nlmixr2(.fastm_one_cmt, d, "foceif", foceiControl(print = 0L, covMethod = "")))
    fRef <- suppressMessages(nlmixr2(.fastm_one_cmt, d, "focei", foceiControl(print = 0L, covMethod = "", fast = TRUE)))
    expect_equal(fF$objf, fRef$objf, tolerance = 0.02)
    expect_equal(unname(fixef(fF)), unname(fixef(fRef)), tolerance = 1e-2)
    # both consumed the analytic gradient and default to lbfgsb3c
    expect_match(fF$extra, "outer: lbfgsb3c; grad: analytic")
    expect_match(fRef$extra, "outer: lbfgsb3c; grad: analytic")
  })

  test_that("all-residual-fixed fast fits use the analytic gradient (omega-only outer)", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    d <- nlmixr2data::theo_sd
    # every residual-error parameter fixed: the mu families profile the
    # structural thetas out too, so the outer problem is omega-only.  The
    # analytic gradient must still engage (it previously fell back to FD, which
    # re-solves the ODE per omega step); rx_r_ is read from the solve, nsg = 0.
    addFix <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
            add.sd <- fix(0.7) })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
              cp <- center / v; cp ~ add(add.sd) })
    }
    combFix <- function() {
      ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
            add.sd <- fix(0.7); prop.sd <- fix(0.1) })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
              d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
              cp <- center / v; cp ~ add(add.sd) + prop(prop.sd) })
    }
    for (.mn in c("addFix", "combFix")) {
      .m <- get(.mn)
      for (.e in c("foceif", "ifoceif", "mfoceif")) {
        .ctl <- switch(.e, foceif = foceiControl, ifoceif = ifoceiControl, mfoceif = mfoceiControl)
        fF <- suppressWarnings(suppressMessages(nlmixr2(.m, d, .e, .ctl(print = 0L, covMethod = ""))))
        ## cached fast=FALSE reference -- see helper-gradref.R
        f0 <- .numRef(paste0("fit-fd-resfix-", .mn, "-", .e), function()
          list(objf = suppressWarnings(suppressMessages(nlmixr2(.m, d, sub("f$", "", .e),
                 .ctl(print = 0L, covMethod = "", fast = FALSE))))$objf))
        expect_match(fF$extra, "grad: analytic", info = .e)
        expect_equal(fF$objf, f0$objf, tolerance = 0.2, info = .e)
      }
    }
  })
})
