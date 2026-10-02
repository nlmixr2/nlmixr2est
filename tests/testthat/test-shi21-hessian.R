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
    .ctl <- nlmControl(
      print = 0L, solveType = "hessian", optimHessType = "central",
      shi21maxHess = 1L, calcTables = FALSE
    )
    .withNlmProblem(.mod, .d, .ctl, function(x) {
      nlmSolveGradHess(x)
      # nlmSolveGradHess() works on R's own vector, which for nlminb is the
      # optimizer's iterate.  The second call used to return with E0 moved to
      # E0 - h.
      x0 <- x + 0
      r <- nlmSolveGradHess(x)
      expect_equal(x, x0, tolerance = 1e-12)
      .h <- attr(r, "hessian")
      # E0's column has no usable difference; the others do
      expect_identical(.h[1, 1], 0)
      expect_true(all(is.finite(.h)))
      expect_true(all(.h[-1, -1] != 0))
    })
  })
})
