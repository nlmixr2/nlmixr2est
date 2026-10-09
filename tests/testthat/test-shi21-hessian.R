nmTest({
  # shi21Hessian() (src/shi21.cpp) builds a Hessian by finite differences of a
  # gradient.  Two callers share it: nlmCalcHessian() (the nlm and nlminb
  # Hessians, est="trust") and calcEtaHessian() (the FOCEi inner Hessian of a
  # non-normal endpoint).

  # The nlm problem for .mod, set up the way a fit sets it up; the scaled
  # starting point is returned and the problem is freed when `f` returns.
  .withNlmProblem <- function(mod, data, control, f) {
    .ui <- rxode2::rxode2(mod)
    .ret <- new.env(parent = emptyenv())
    .foceiPreProcessData(data, .ret, .ui, control$rxControl)
    .p <- setNames(.ui$nlmParIni, .ui$nlmParName)
    on.exit(.nlmFreeEnv())
    .env <- .nlmSetupEnv(.p, .ui, .ret$dataSav, .ui$nlmSensModel, control)
    f(.env$par.ini + 0)
  }

  test_that("a Hessian column with both legs non-finite leaves theta where it was", {
    skip_on_cran()
    # exp() overflows to Inf once |E0 - 0.5| > 1.31e-3, so the gradient in E0
    # is infinite on both sides of E0 = 0.5 at any step bigger than that.
    .mod <- function() {
      ini({
        E0 <- 0.5
        Em <- 0.5
        E50 <- 2
      })
      model({
        v <- E0 + Em * time / (E50 + time)
        ll(bin) ~ DV * v - log(1 + exp(v)) - exp(1e9 * ((E0 - 0.5)^2 - 1e-6))
      })
    }
    .d <- data.frame(ID = 1L, TIME = seq(0.1, 10, length.out = 60), AMT = 0, EVID = 0L)
    .d$DV <- as.integer(seq_len(nrow(.d)) %% 3 == 0)
    # shi21maxHess = 1: the first call's one-iteration step search finds no
    # finite E0 leg and keeps its starting step (2.6e-2 on the scaled theta),
    # so the second call differences E0 at that step, outside the window.
    .ctl <- nlmControl(print = 0L, solveType = "hessian", optimHessType = "central", shi21maxHess = 1L)
    .withNlmProblem(.mod, .d, .ctl, function(x) {
      nlmSolveGradHess(x)
      # nlmSolveGradHess() works on R's own vector, which for nlminb is the
      # optimizer's iterate: the Hessian has to leave it exactly as it was.
      x0 <- x + 0
      r <- nlmSolveGradHess(x)
      expect_identical(x, x0)
      .h <- attr(r, "hessian")
      # E0's column has no usable difference; the others do
      expect_identical(.h[1, 1], 0)
      expect_true(all(is.finite(.h)))
      expect_true(all(.h[-1, -1] != 0))
    })
  })

  test_that("the Hessian matches the gradient's Jacobian, searched and cached", {
    skip_on_cran()
    # Cubic in a and quadratic otherwise: the gradient is quadratic, so a
    # central difference of it is exact at any step, while a forward one is off
    # in the (a, a) cell by an amount set by a's step.
    .mod <- function() {
      ini({
        a <- 0.3
        b <- -0.2
        c <- 0.7
      })
      model({
        v <- a + b * time
        ll(bin) ~ DV * v + 0.5 * a^3 - 2 * (a - 0.5)^2 - a * b - 3 * (b + 1)^2 + b * c + 0.25 * a * c - 0.5 * c^2
      })
    }
    .d <- data.frame(ID = 1L, TIME = seq(0.1, 10, length.out = 20), AMT = 0, EVID = 0L)
    .d$DV <- as.integer(seq_len(nrow(.d)) %% 2 == 0)
    .hess <- list()
    .oracle <- NULL
    for (.type in c("central", "forward")) {
      .ctl <- nlmControl(print = 0L, solveType = "hessian", optimHessType = .type)
      .withNlmProblem(.mod, .d, .ctl, function(x) {
        .gr <- function(p) attr(nlmSolveGradR(p), "gradient")
        .h <- 1e-3
        .oracle <<- vapply(
          seq_along(x),
          function(k) {
            .e <- replace(numeric(length(x)), k, .h)
            (.gr(x + .e) - .gr(x - .e)) / (2 * .h)
          },
          numeric(length(x))
        )
        # first call searches the steps, second reuses them
        .hess[[.type]] <<- lapply(1:2, function(.i) attr(nlmSolveGradHess(x + 0), "hessian"))
      })
    }
    # cross terms are nonzero, so a swapped column or row would show
    expect_true(all(.oracle[upper.tri(.oracle)] != 0))
    for (.i in 1:2) {
      expect_equal(.hess$central[[.i]], .oracle, tolerance = 1e-6, info = .i)
    }
    # forward differs from the oracle in exactly the (a, a) cell
    .off <- abs(.hess$forward[[1]] - .oracle) > 1e-6 * abs(.oracle)
    expect_equal(sum(.off), 1L)
    expect_true(.off[1, 1])
    # the step cached for the second call is the one the search settled on
    expect_equal(.hess$forward[[2]], .hess$forward[[1]], tolerance = 1e-10)
  })

  test_that("a gradient component the column does not touch leaves its step alone (#1188)", {
    skip_on_cran()
    # g_c does not depend on a, so its third difference along a is roundoff; the
    # legacy ratio let it pin the search ratio near 0 and cap a's step at hMax.
    .mod <- function() {
      ini({
        a <- 0.3
        b <- -0.2
        c <- 0.7
      })
      model({
        v <- a + b * time
        ll(bin) ~ DV * v - log(1 + exp(v)) - exp(4 * a) - 0.1 * exp(0.5 * b) - 0.5 * c^2 + 0.01 * c^4
      })
    }
    .d <- data.frame(ID = 1L, TIME = seq(0.1, 10, length.out = 20), AMT = 0, EVID = 0L)
    .d$DV <- as.integer(seq_len(nrow(.d)) %% 2 == 0)
    .ctl <- nlmControl(print = 0L, solveType = "hessian", optimHessType = "central")
    .old <- .shi21RatioCensor()
    on.exit(.shi21RatioCensor(.old))
    # .nlmSetupEnv() reads the option
    .err <- vapply(
      c("legacy", "detected"),
      function(.type) {
        withr::local_options(list(nlmixr2est.shi21RatioCensor = .type))
        .withNlmProblem(.mod, .d, .ctl, function(x) {
          .gr <- function(p) attr(nlmSolveGradR(p + 0), "gradient")
          .oracle <- numDeriv::jacobian(.gr, x + 0)
          .oracle <- (.oracle + t(.oracle)) / 2
          .h <- attr(nlmSolveGradHess(x + 0), "hessian")
          # cross terms can be exactly 0, so scale by the diagonal
          max(abs(.h - .oracle) / sqrt(abs(outer(diag(.oracle), diag(.oracle)))))
        })
      },
      numeric(1)
    )
    expect_gt(.err[["legacy"]], 1)
    expect_lt(.err[["detected"]], 1e-3)

    # a fit reads the option: nlminb on the legacy Hessian stops well short
    .fit <- function(type) {
      withr::with_options(list(nlmixr2est.shi21RatioCensor = type), {
        suppressMessages(nlmixr2(.mod, .d, est = "nlminb", control = nlminbControl(print = 0L)))$objf
      })
    }
    .ofv <- vapply(c("legacy", "detected"), .fit, numeric(1))
    expect_gt(.ofv[["legacy"]] - .ofv[["detected"]], 1)
  })

  test_that("shi21HessRefresh re-searches a step once theta leaves its search span (#1175)", {
    skip_on_cran()
    # different curvature along each coordinate, so the searched steps differ
    .mod <- function() {
      ini({
        a <- 0.3
        b <- -0.2
        c <- 0.7
      })
      model({
        v <- a + b * time
        ll(bin) ~ DV * v - log(1 + exp(v)) - exp(4 * a) - 0.1 * exp(0.5 * b) - 0.5 * c^2 + 0.01 * c^4
      })
    }
    .d <- data.frame(ID = 1L, TIME = seq(0.1, 10, length.out = 20), AMT = 0, EVID = 0L)
    .d$DV <- as.integer(seq_len(nrow(.d)) %% 2 == 0)
    for (.type in c("central", "forward")) {
      .ctl <- nlmControl(
        print = 0L,
        solveType = "hessian",
        optimHessType = .type,
        shi21HessRefresh = TRUE
      )
      .withNlmProblem(.mod, .d, .ctl, function(x) {
        nlmSolveGradHess(x + 0)
        .i0 <- .nlmHessStepInfo()
        .h <- .i0$step
        expect_equal(.i0$nSearch, 3L)
        expect_true(all(.h > 0))
        # the steps differ, so a gate mixing one coordinate's move with another's
        # step would show below
        expect_gt(max(.h) / min(.h), 1.2)
        .span <- if (.type == "central") 3 else 4
        # every coordinate inside its own step's span: steps reused
        nlmSolveGradHess(x + 0)
        nlmSolveGradHess(x + 0.97 * .span * .h)
        expect_identical(.nlmHessStepInfo()[c("step", "nSearch")], .i0[c("step", "nSearch")], info = .type)
        # only the smallest-step coordinate leaves its span: every step re-searched
        .k <- which.min(.h)
        nlmSolveGradHess(replace(x, .k, x[.k] + 1.03 * .span * .h[.k]))
        expect_equal(.nlmHessStepInfo()$nSearch, 6L, info = .type)
        # every coordinate past its span
        nlmSolveGradHess(x + 1.1 * .span * .h)
        expect_equal(.nlmHessStepInfo()$nSearch, 9L, info = .type)
      })
    }
  })

  test_that("a re-searched Hessian step starts from the old step (#1175)", {
    skip_on_cran()
    .mod <- function() {
      ini({
        a <- 0.3
        b <- -0.2
        c <- 0.7
      })
      model({
        v <- a + b * time
        ll(bin) ~ DV * v - log(1 + exp(v)) - exp(4 * a) - 0.1 * exp(0.5 * b) - 0.5 * c^2 + 0.01 * c^4
      })
    }
    .d <- data.frame(ID = 1L, TIME = seq(0.1, 10, length.out = 20), AMT = 0, EVID = 0L)
    .d$DV <- as.integer(seq_len(nrow(.d)) %% 2 == 0)
    # two search iterations: a step that keeps growing ends one growth past where
    # its search started, so a warm re-search ends past the old step while a
    # cold one would land on it again.  The legacy ratio keeps every step growing
    # here (#1188).
    .old <- .shi21RatioCensor()
    on.exit(.shi21RatioCensor(.old))
    withr::local_options(list(nlmixr2est.shi21RatioCensor = "legacy"))
    .ctl <- nlmControl(
      print = 0L,
      solveType = "hessian",
      optimHessType = "central",
      shi21maxHess = 2L,
      shi21HessRefresh = TRUE
    )
    .withNlmProblem(.mod, .d, .ctl, function(x) {
      nlmSolveGradHess(x + 0)
      .h0 <- .nlmHessStepInfo()$step
      nlmSolveGradHess(x + 1.1 * 3 * .h0)
      .i1 <- .nlmHessStepInfo()
      expect_equal(.i1$nSearch, 6L)
      expect_true(all(.i1$step > .h0))
    })
  })

  test_that("nlm and nlminb re-search the Hessian steps only under shi21HessRefresh (#1175)", {
    skip_on_cran()
    .mod <- function() {
      ini({
        a <- 0.3
        b <- -0.2
      })
      model({
        v <- a + b * time
        ll(bin) ~ DV * v - log(1 + exp(v))
      })
    }
    .d <- data.frame(ID = 1L, TIME = seq(0.1, 10, length.out = 40), AMT = 0, EVID = 0L)
    .d$DV <- as.integer(seq_len(nrow(.d)) %% 3 != 0 & seq_len(nrow(.d)) < 30)
    for (.refresh in c(TRUE, FALSE)) {
      .ctl <- nlmControl(
        print = 0L,
        solveType = "hessian",
        optimHessType = "central",
        shi21HessRefresh = .refresh
      )
      # without the refresh, the two steps are searched once, at the start
      .check <- if (.refresh) {
        function() expect_gt(.nlmHessStepInfo()$nSearch, 2L)
      } else {
        function() expect_equal(.nlmHessStepInfo()$nSearch, 2L)
      }
      .withNlmProblem(.mod, .d, .ctl, function(x) {
        stats::nlm(function(p) nlmSolveGradHess(p), x + 0, check.analyticals = FALSE)
        .check()
      })
      .withNlmProblem(.mod, .d, .ctl, function(x) {
        stats::nlminb(
          x + 0,
          function(p) nlminbFunC(p, 1L),
          gradient = function(p) nlminbFunC(p, 2L),
          hessian = function(p) nlminbFunC(p, 3L)
        )
        .check()
      })
    }
  })

  test_that("optimHessType='richardson' extrapolates the central Hessian (#1175)", {
    skip_on_cran()
    skip_if_not_installed("numDeriv")
    .mod <- function() {
      ini({
        tka <- log(0.5)
        tcl <- log(5)
        tv <- log(50)
        add.sd <- 2
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        cp <- linCmt()
        cp ~ add(add.sd)
      })
    }
    .hess <- list()
    .oracle <- NULL
    for (.type in c("central", "richardson")) {
      .ctl <- nlmControl(print = 0L, solveType = "hessian", optimHessType = .type)
      .withNlmProblem(.mod, nlmixr2data::theo_sd, .ctl, function(x) {
        .gr <- function(p) attr(nlmSolveGradR(p), "gradient")
        .oracle <<- numDeriv::jacobian(.gr, x + 0)
        .oracle <<- 0.5 * (.oracle + t(.oracle))
        nlmSolveGradHess(x + 0)
        .i0 <- .nlmHessStepInfo()
        # the second call reuses the searched steps
        .hess[[.type]] <<- attr(nlmSolveGradHess(x + 0), "hessian")
        .i1 <- .nlmHessStepInfo()
        expect_equal(.i1$nSearch, 4L)
        # four gradient solves per coordinate for Richardson, two for central
        expect_equal(.i1$nGrad - .i0$nGrad, if (.type == "central") 8L else 16L, info = .type)
      })
    }
    .err <- vapply(.hess, function(h) max(abs(h - .oracle)) / max(abs(.oracle)), numeric(1))
    expect_lt(.err[["richardson"]], 1e-6)
    expect_lt(.err[["richardson"]], .err[["central"]] / 100)
  })

  test_that("a non-normal-endpoint FOCEi fit reports llikObs at its final ETAs", {
    skip_on_cran()
    # A dnorm() endpoint sets needOptimHess: the inner Hessian is a finite
    # difference of the eta gradient, and every leg re-solves the subject,
    # rewriting its per-observation log-likelihoods.  The fit still has to
    # report them at its ETAs, not at the last leg.
    .one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- log(2.7)
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        cp <- linCmt()
        cp ~ add(add.sd) + dnorm()
      })
    }
    .check <- function(control) {
      .fit <- suppressMessages(suppressWarnings(
        nlmixr2(.one.cmt, nlmixr2data::theo_sd, "focei", control = control)
      ))
      .ll <- .fit$llikObs
      .ll <- .ll[!is.na(.ll)] # dose records
      # the table's IPRED is solved at the fit's ETAs
      expect_equal(
        .ll,
        dnorm(.fit$DV, .fit$IPRED, .fit$theta[["add.sd"]], log = TRUE),
        tolerance = 1e-10
      )
    }
    # central and forward differences (the final objective re-searches the
    # steps), and the trust inner optimizer
    .check(foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L))
    .check(foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L, optimHessCovType = "forward"))
    .check(foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L, innerOpt = "trust", hessianMethod = "fd"))
  })
})
