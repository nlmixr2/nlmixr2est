nmTest({
  # foceiControl(etaPolish=): the analytic gradient is exact only at the eta mode, but the
  # inner optimizer stops at its own tolerance.  "none" leaves the etas there; "newton1"
  # takes one exact-Hessian Newton step; "newton" iterates to round-off, so the objective
  # and the gradient are taken at the same eta.
  .ebePolishModel <- function() {
    ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; wtcl <- 0.4
      eta.cl ~ 0.08; eta.v ~ 0.02; add.sd <- 0.7; prop.sd <- 0.08 })
    model({
      ka <- exp(tka); cl <- exp(tcl + eta.cl) * (WT / 70)^wtcl; v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v; cp ~ add(add.sd) + prop(prop.sd)
    })
  }
  .ebePolishFit <- function(est, mode, tight = FALSE, ...) {
    ctl <- list(
      print = 0L,
      calcTables = FALSE,
      covMethod = "",
      fast = TRUE,
      maxOuterIterations = 0L,
      rxControl = rxode2::rxControl(atol = 1e-11, rtol = 1e-11),
      ...
    )
    if (!is.null(mode)) {
      ctl$etaPolish <- mode
    }
    if (tight) {
      ctl <- c(
        ctl,
        list(innerOpt = "trust", trustPolish = TRUE, trustFterm = 1e-15, trustMterm = 1e-15, maxInnerIterations = 2000L)
      )
    }
    suppressMessages(suppressWarnings(nlmixr2(.ebePolishModel, nlmixr2data::theo_sd, est, do.call(foceiControl, ctl))))
  }
  # gradient of the reported objective by Richardson extrapolation of central differences
  .ebePolishFd <- function(fit, nm) {
    th <- fit$theta
    cd <- function(n, h) {
      hh <- h * max(abs(th[[n]]), 1)
      up <- dn <- th[n]
      up[] <- th[[n]] + hh
      dn[] <- th[[n]] - hh
      (nlmixr2ObjGrad(fit, theta = up)$objective - nlmixr2ObjGrad(fit, theta = dn)$objective) / (2 * hh)
    }
    vapply(nm, function(n) (4 * cd(n, 5e-4) - cd(n, 1e-3)) / 3, numeric(1))
  }
  .ebePolishNm <- c("tcl", "tv", "add.sd")

  test_that("etaPolish control: default, validation, option", {
    expect_identical(foceiControl()$etaPolish, "none")
    expect_identical(foceiControl(etaPolish = "newton1")$etaPolish, "newton1")
    expect_identical(agqControl(etaPolish = "newton")$etaPolish, "newton")
    expect_identical(laplaceControl(etaPolish = "newton1")$etaPolish, "newton1")
    expect_error(foceiControl(etaPolish = "chord"))
    withr::with_options(list(nlmixr2est.etaPolish = "newton1"), expect_identical(foceiControl()$etaPolish, "newton1"))
  })

  .cases <- list(list(est = "focei", args = list()), list(est = "foce", args = list(foce = "nonmem")))
  for (.case in .cases) {
    test_that(paste("etaPolish", .case$est, ": default is none; newton1 and newton polish the eta mode"), {
      skip_on_cran()
      skip_if_not_installed("nlmixr2data")
      fit <- function(mode, tight = FALSE) do.call(.ebePolishFit, c(list(.case$est, mode, tight), .case$args))
      def <- nlmixr2ObjGrad(fit(NULL))
      non <- nlmixr2ObjGrad(fit("none"))
      # "none" is the default, bit for bit
      expect_identical(non$objective, def$objective)
      expect_identical(non$gradient, def$gradient)
      ref <- nlmixr2ObjGrad(fit("newton", TRUE)) # tight inner: the exact mode
      nw <- nlmixr2ObjGrad(fit("newton"))
      nw1 <- nlmixr2ObjGrad(fit("newton1"))
      dObj <- function(r) abs(r$objective - ref$objective)
      dGr <- function(r) max(abs(r$gradient - ref$gradient) / pmax(abs(ref$gradient), 1))
      # measured on theo (relative objective / gradient, default vs tight inner): none
      # 2e-7 / 2e-6 (FOCE 9e-7 / 1e-4), newton1 2e-11 / 5e-11 (3 etas: 4e-9 / 5e-9), newton
      # 4e-14 / 2e-13
      expect_gt(dObj(non), 1e-7)
      expect_lt(dObj(nw), 1e-8)
      expect_lt(dGr(nw), 1e-8)
      expect_lt(dObj(nw1), 1e-5)
      expect_lt(dGr(nw1), 1e-5)
      expect_lt(dObj(nw1), dObj(non))
      expect_lt(dGr(nw1), dGr(non))
    })
  }

  test_that("etaPolish: the analytic gradient is the derivative of the reported objective", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    nm <- c("tcl", "add.sd")
    err <- function(mode) {
      f <- .ebePolishFit("focei", mode)
      g <- nlmixr2ObjGrad(f)$gradient[nm]
      fd <- .ebePolishFd(f, nm)
      max(abs(g - fd)) / max(abs(fd))
    }
    eNone <- err("none")
    eOne <- err("newton1")
    eNew <- err("newton")
    # measured: none 4e-4, newton1 1e-9, newton 1e-9 (the finite difference's own noise)
    expect_gt(eNone, 1e-5)
    expect_lt(eOne, 1e-6)
    expect_lt(eNew, 1e-6)
  })
})
