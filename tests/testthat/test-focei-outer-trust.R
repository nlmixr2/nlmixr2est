test_that("outerOpt='trust' validates its own control", {
  expect_error(foceiControl(outerOpt = "trust", outerTrustHessian = "analytic"),
               "requires fast")
  expect_no_error(foceiControl(outerOpt = "trust", fast = TRUE,
                               outerTrustHessian = "analytic"))
  expect_error(foceiControl(outerTrustRinit = 0), "must be > 0")
  expect_error(foceiControl(outerTrustRmax = 0), "must be > 0")
  expect_error(foceiControl(outerTrustRinit = 2, outerTrustRmax = 1),
               "cannot be larger")
  expect_error(foceiControl(outerTrustRelStep = 0), "must be > 0")
  expect_error(foceiControl(outerTrustRestarts = -1L))
  # a derivative-based outer optimizer must NOT trip the derivative-free
  # fast= downgrade that bobyqa/uobyqa/newuoa do
  expect_true(foceiControl(outerOpt = "trust", fast = TRUE)$fast)
  expect_equal(foceiControl(outerOpt = "trust")$outerOptTxt, "trust")
})

test_that("the Newton decrement gate reads a trust result", {
  # positive definite Hessian: 0.5 * g' H^-1 g
  expect_equal(.trustOuterDecrement(list(gradient = c(1, 2),
                                         hessian = diag(c(2, 8)))),
               0.5 * (1 / 2 + 4 / 8))
  # an indefinite Hessian at a reported minimum is not a minimum
  expect_true(is.na(.trustOuterDecrement(list(gradient = c(1, 0),
                                              hessian = diag(c(-1, 1))))))
  expect_true(is.na(.trustOuterDecrement(list(gradient = c(NA_real_, 0),
                                              hessian = diag(2)))))
  expect_true(is.na(.trustOuterDecrement(list(gradient = c(1, 0),
                                              hessian = NULL))))
})

test_that("the curvature supplier falls back when the analytic Hessian declines", {
  .fn <- function(x) sum(x^2)
  .gr <- function(x) c(2 * x[1], 2 * x[2])
  .box <- c(-Inf, -Inf)
  .hi <- c(Inf, Inf)
  # fast=FALSE: the analytic Hessian is not available at all.  foceiControl()
  # refuses the pinned request, but `fast` can still be downgraded after that, so
  # the runtime demotes with a warning rather than aborting the fit.
  expect_warning(
    .c <- .trustOuterCurvature(
      list(outerTrustHessian = "analytic", fast = FALSE,
           hessian = function(x, relStep) diag(2)),
      .fn, .gr, 1e-3, .box, .hi
    ),
    "needs fast"
  )
  expect_equal(.c$hessian(c(1, 1), .gr(c(1, 1))), diag(2))
  expect_equal(.c$calls, 0L)
  expect_false(.c$fallback)

  # available, then refused mid-run: one warning, and every later call goes to
  # BFGS rather than paying the failed probe again
  .n <- 0L
  .ctl <- list(outerTrustHessian = "analytic", fast = TRUE,
               hessian = function(x, relStep) {
                 .n <<- .n + 1L
                 stop("analytical outer Hessian unavailable (status -4)")
               })
  .c <- .trustOuterCurvature(.ctl, .fn, .gr, 1e-3, .box, .hi)
  expect_warning(.c$hessian(c(1, 1), .gr(c(1, 1))), "continues with BFGS")
  expect_true(.c$fallback)
  expect_equal(.c$calls, 1L)
  expect_silent(.c$hessian(c(1.1, 1), .gr(c(1.1, 1))))
  expect_equal(.n, 1L)
  expect_equal(.c$calls, 1L)

  # "fd" declines when neither difference direction fits in the box; the
  # supplier answers from BFGS instead of returning NULL to the optimizer
  .c <- .trustOuterCurvature(
    list(outerTrustHessian = "fd", fast = FALSE), .fn, .gr, 1e-3, c(1, 1), c(1, 1)
  )
  expect_equal(.c$hessian(c(1, 1), .gr(c(1, 1))), diag(2))

  # outerTrustRelStep must actually reach the analytic entry -- foceiOuterH()
  # takes it as `relStep` and silently keeps its own default otherwise
  .seen <- NULL
  .c <- .trustOuterCurvature(
    list(outerTrustHessian = "analytic", fast = TRUE,
         hessian = function(x, relStep) { .seen <<- relStep; diag(2) }),
    .fn, .gr, 5e-4, .box, .hi
  )
  .c$hessian(c(1, 1), .gr(c(1, 1)))
  expect_equal(.seen, 5e-4)
})

test_that("the finite-difference curvature settles every point it reads", {
  # the gradient callback warm-starts from the last evaluation, so a probe read
  # without settling it first returns the gradient at a stale conditional mode
  .seen <- character()
  .fn <- function(x) {
    .seen <<- c(.seen, paste0("fn:", paste(signif(x, 8), collapse = ",")))
    sum(x^2)
  }
  .gr <- function(x) {
    .seen <<- c(.seen, paste0("gr:", paste(signif(x, 8), collapse = ",")))
    c(2 * x[1], 2 * x[2])
  }
  .fd <- .trustOuterFd(.fn, .gr, 1e-3, c(-Inf, -Inf), c(Inf, Inf))
  .g0 <- .gr(c(1, 1))
  .seen <- character()
  .h <- .fd(c(1, 1), .g0)
  expect_equal(.h, diag(c(2, 2)), tolerance = 1e-6)
  # every gr the difference took was preceded by an fn at the SAME point
  .at <- which(startsWith(.seen, "gr:"))
  expect_length(.at, 2L)
  for (.i in .at) {
    expect_gt(.i, 1L)
    expect_identical(sub("^fn:", "", .seen[.i - 1L]), sub("^gr:", "", .seen[.i]))
  }
  # ... and the point is settled again on the way out
  expect_identical(.seen[length(.seen)], "fn:1,1")
})

test_that("outerOpt='trust' fits and consumes the analytic outer Hessian", {
  skip_on_cran()
  model <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
          eta.cl ~ 0.3; add.sd <- 0.7 })
    model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
            d/dt(depot) <- -ka * depot
            d/dt(center) <- ka * depot - cl / v * center
            cp <- center / v
            cp ~ add(add.sd) })
  }
  d <- nlmixr2data::theo_sd
  ctl <- function(...) {
    foceiControl(print = 0L, calcTables = FALSE, covMethod = "",
                 outerOpt = "trust", ...)
  }
  fitA <- .nlmixr(model, d, "focei", ctl(fast = TRUE))
  # The counter is what proves the analytic Hessian ran; equal objectives alone
  # cannot tell it from the quasi-Newton fallback.
  expect_gt(fitA$env$optReturn$hessianEvaluations, 0L)
  expect_false(fitA$env$optReturn$hessianFallback)
  expect_true(is.finite(fitA$objf))

  fitB <- .nlmixr(model, d, "focei", ctl(fast = TRUE, outerTrustHessian = "bfgs"))
  expect_equal(fitB$env$optReturn$hessianEvaluations, 0L)
  expect_equal(fitB$objf, fitA$objf, tolerance = 1e-3)

  # ... and it reaches the same optimum as the shipping outer optimizer
  fitN <- .nlmixr(model, d, "focei",
                  foceiControl(print = 0L, calcTables = FALSE, covMethod = "",
                               outerOpt = "nlminb", fast = TRUE))
  expect_equal(fitA$objf, fitN$objf, tolerance = 1e-3)
  expect_equal(unname(fixef(fitA)), unname(fixef(fitN)), tolerance = 1e-2)
})
