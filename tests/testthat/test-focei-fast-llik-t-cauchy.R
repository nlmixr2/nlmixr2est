nmTest({
  # fast = TRUE built the outer model and innerHess2 with rx_pred_ = f for a t()/cauchy()
  # endpoint, while the ll() machinery reads rx_pred_ as the log-density: the objective
  # was wrong (and the gradient with it).  The fast objective must equal the fast = FALSE one.
  .m <- function(err, sd = quote(add.sd <- 0.7)) {
    f <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        .sd
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        .err
      })
    }
    body(f) <- eval(substitute(substitute(b, list(.err = err, .sd = sd)), list(b = body(f))))
    f
  }
  .fit <- function(m, est, fast, ...) {
    ctl <- if (est == "laplace") laplaceControl else foceiControl
    suppressMessages(suppressWarnings(nlmixr2(
      m,
      nlmixr2data::theo_sd,
      est,
      ctl(print = 0L, covMethod = "", calcTables = FALSE, fast = fast, ...)
    )))
  }
  .errs <- list(
    t = quote(cp ~ add(add.sd) + dt(5)),
    cauchy = quote(cp ~ add(add.sd) + dcauchy())
  )
  for (.n in names(.errs)) {
    for (.est in c("focei", "laplace")) {
      test_that(paste0(.n, " endpoint: fast objective equals the non-fast one (", .est, ")"), {
        skip_on_cran()
        skip_if_not_installed("nlmixr2data")
        .m1 <- .m(.errs[[.n]])
        # fast = TRUE uses the exact inner Hessian at a tightly converged eta-hat.  The
        # fast = FALSE eta-hat stops at sigdig = 3 (objective low by ~2e-3 for t and ~7e-3
        # for cauchy, an eta-hat convergence gap, not a Hessian one), so compare at
        # sigdig = 8, where the two agree to ~1e-6; the same value comes from an
        # independent Laplace computation (FD Hessian, optim-converged eta).
        .a <- .fit(.m1, .est, FALSE, maxOuterIterations = 0L, sigdig = 8)
        .b <- .fit(.m1, .est, TRUE, maxOuterIterations = 0L, sigdig = 8)
        expect_equal(.b$objf, .a$objf, tolerance = 1e-6)
      })
    }
  }
  # an analytic gradient that declines (nAnalyticGradDirect == 0) must fail the test
  .gradCheck <- function(err, extra = NULL, ...) {
    .m1 <- if (is.null(extra)) .m(err) else .m(err, extra)
    .f <- .fit(.m1, "focei", TRUE, maxOuterIterations = 0L, ...)
    .th <- .f$theta[c("tka", "tcl", "tv", "add.sd")]
    .th["tcl"] <- .th["tcl"] + 0.05
    .th["add.sd"] <- .th["add.sd"] * 1.1
    .om <- .f$omega * 1.1
    .g <- suppressWarnings(nlmixr2ObjGrad(.f, theta = .th, omega = .om))
    expect_gt(attr(.g, "nAnalyticGradDirect"), 0L)
    .obj <- function(th) suppressWarnings(nlmixr2ObjGrad(.f, theta = th, omega = .om))$objective
    for (.p in names(.th)) {
      .h <- 1e-3
      .tp <- .th
      .tm <- .th
      .tp[.p] <- .tp[.p] + .h
      .tm[.p] <- .tm[.p] - .h
      expect_equal(.g$gradient[[.p]], (.obj(.tp) - .obj(.tm)) / (2 * .h), tolerance = 5e-3)
    }
  }
  test_that("t endpoint: analytic gradient (incl. sigma) matches central differences", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .gradCheck(.errs$t)
  })
  # a fixed transform: the closed form carries it, so the outer model must not announce
  # one (rx_predf_ is a log-density); sigdig = 8 keeps the central-difference noise down
  test_that("t endpoint with a fixed boxCox: analytic gradient matches central differences", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .gradCheck(
      quote(cp ~ add(add.sd) + boxCox(lam) + dt(5)),
      quote({
      add.sd <- 0.7
      lam <- fix(0.5)
    }),
      sigdig = 8
    )
  })
})
