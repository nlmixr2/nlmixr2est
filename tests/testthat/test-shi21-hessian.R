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
    # Quadratic in the thetas, so the Hessian is constant and any correct
    # difference of the gradient recovers it to rounding.
    .mod <- function() {
      ini({
        a <- 0.3
        b <- -0.2
        c <- 0.7
      })
      model({
        v <- a + b * time
        ll(bin) ~ DV * v - 2 * (a - 0.5)^2 - a * b - 3 * (b + 1)^2 + b * c + 0.25 * a * c - 0.5 * c^2
      })
    }
    .d <- data.frame(ID = 1L, TIME = seq(0.1, 10, length.out = 20), AMT = 0, EVID = 0L)
    .d$DV <- as.integer(seq_len(nrow(.d)) %% 2 == 0)
    for (.type in c("central", "forward")) {
      .ctl <- nlmControl(print = 0L, solveType = "hessian", optimHessType = .type)
      .withNlmProblem(.mod, .d, .ctl, function(x) {
        .gr <- function(p) attr(nlmSolveGradR(p), "gradient")
        .h <- 1e-3
        .oracle <- vapply(
          seq_along(x),
          function(k) {
            .e <- replace(numeric(length(x)), k, .h)
            (.gr(x + .e) - .gr(x - .e)) / (2 * .h)
          },
          numeric(length(x))
        )
        # cross terms are nonzero, so a swapped column or row would show
        expect_true(all(.oracle[upper.tri(.oracle)] != 0))
        # first call searches the steps, second reuses them
        for (.i in 1:2) {
          .hess <- attr(nlmSolveGradHess(x + 0), "hessian")
          expect_equal(.hess, .oracle, tolerance = 1e-6, info = paste(.type, .i))
        }
      })
    }
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
