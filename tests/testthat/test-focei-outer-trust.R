test_that("outerOpt='trust' validates its own control", {
  expect_error(foceiControl(outerOpt = "trust", outerTrustHessian = "analytic"), "requires fast")
  expect_no_error(foceiControl(outerOpt = "trust", fast = TRUE, outerTrustHessian = "analytic"))
  expect_error(foceiControl(outerTrustRinit = 0), "must be > 0")
  expect_error(foceiControl(outerTrustRmax = 0), "must be > 0")
  expect_error(foceiControl(outerTrustRinit = 2, outerTrustRmax = 1), "cannot be larger")
  expect_error(foceiControl(outerTrustRelStep = 0), "must be > 0")
  expect_error(foceiControl(outerTrustRestarts = -1L))
  # a derivative-based outer optimizer must NOT trip the derivative-free
  # fast= downgrade that bobyqa/uobyqa/newuoa do
  expect_true(foceiControl(outerOpt = "trust", fast = TRUE)$fast)
  expect_equal(foceiControl(outerOpt = "trust")$outerOptTxt, "trust")
})

test_that("the trust outer controls survive a control round trip", {
  # A built control is rebuilt from itself (posthoc re-validation, the *f
  # wrappers, rxUiDeparse), so every new argument has to come back as the type
  # and value it went in as -- outerOpt in particular, which is stored as -1L
  # with the function in outerOptFun and recovered by name from outerOptTxt.
  .c <- foceiControl(
    outerOpt = "trust",
    fast = TRUE,
    outerTrustHessian = "fd",
    outerTrustRinit = 0.3,
    outerTrustRmax = 2,
    outerTrustFterm = 1e-7,
    outerTrustMterm = 1e-8,
    outerTrustRelStep = 5e-4,
    outerTrustRestarts = 5L
  )
  .r <- do.call(foceiControl, unclass(.c))
  for (.n in c(
    "outerTrustHessian",
    "outerTrustRinit",
    "outerTrustRmax",
    "outerTrustFterm",
    "outerTrustMterm",
    "outerTrustRelStep",
    "outerTrustRestarts",
    "outerOptTxt"
  )) {
    expect_identical(.r[[.n]], .c[[.n]], info = .n)
  }
  expect_true(is.function(.r$outerOptFun))
  expect_match(paste(deparse(rxode2::rxUiDeparse(.c, "ctl")), collapse = " "), "outerOpt = \"trust\"")
})

test_that("the Newton decrement gate reads a trust result", {
  # positive definite Hessian: 0.5 * g' H^-1 g
  expect_equal(.trustOuterDecrement(list(gradient = c(1, 2), hessian = diag(c(2, 8)))), 0.5 * (1 / 2 + 4 / 8))
  # an indefinite Hessian at a reported minimum is not a minimum
  expect_true(is.na(.trustOuterDecrement(list(gradient = c(1, 0), hessian = diag(c(-1, 1))))))
  expect_true(is.na(.trustOuterDecrement(list(gradient = c(NA_real_, 0), hessian = diag(2)))))
  expect_true(is.na(.trustOuterDecrement(list(gradient = c(1, 0), hessian = NULL))))
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
      list(outerTrustHessian = "analytic", fast = FALSE, hessian = function(x, relStep) diag(2)),
      .fn,
      .gr,
      1e-3,
      .box,
      .hi
    ),
    "needs fast"
  )
  expect_equal(.c$hessian(c(1, 1), .gr(c(1, 1))), diag(2))
  expect_equal(.c$calls, 0L)
  expect_false(.c$fallback)

  # available, then refused mid-run: one warning, and every later call goes to
  # BFGS rather than paying the failed probe again
  .n <- 0L
  .ctl <- list(outerTrustHessian = "analytic", fast = TRUE, hessian = function(x, relStep) {
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
    list(outerTrustHessian = "fd", fast = FALSE),
    .fn,
    .gr,
    1e-3,
    c(1, 1),
    c(1, 1)
  )
  expect_equal(.c$hessian(c(1, 1), .gr(c(1, 1))), diag(2))

  # outerTrustRelStep must actually reach the analytic entry -- foceiOuterH()
  # takes it as `relStep` and silently keeps its own default otherwise
  .seen <- NULL
  .c <- .trustOuterCurvature(
    list(outerTrustHessian = "analytic", fast = TRUE, hessian = function(x, relStep) {
      .seen <<- relStep
      diag(2)
    }),
    .fn,
    .gr,
    5e-4,
    .box,
    .hi
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

test_that("the trust driver hands its control through to the region and curvature", {
  # The unit tests above exercise the helpers directly, so they cannot see the
  # driver reading the wrong control field.  Drive .trustOuter() itself over a
  # quadratic whose answer is known and watch what reaches each seam.
  .seen <- list()
  .control <- list(
    fast = TRUE,
    sigdig = 3,
    maxOuterIterations = 50L,
    outerTrustHessian = "analytic",
    outerTrustRelStep = 5e-4,
    outerTrustRinit = 0.4,
    outerTrustRmax = 3.2,
    outerTrustRestarts = 2L,
    outerTrustFterm = 1e-11,
    outerTrustMterm = 1e-11,
    hessian = function(x, relStep) {
      .seen$relStep <<- relStep
      diag(c(2, 8))
    }
  )
  .ret <- .trustOuter(
    c(3, -2),
    fn = function(x) x[1]^2 + 4 * x[2]^2,
    gr = function(x) c(2 * x[1], 8 * x[2]),
    lower = c(-Inf, -Inf),
    upper = c(Inf, Inf),
    control = .control
  )
  expect_equal(.ret$x, c(0, 0), tolerance = 1e-6)
  expect_equal(.ret$convergence, 0L)
  expect_equal(.seen$relStep, 5e-4)
  expect_gt(.ret$hessianEvaluations, 0L)
  expect_false(.ret$hessianFallback)
  # a true Newton step from a quadratic lands in one iteration, so the restart
  # budget is untouched and the decrement gate is satisfied
  expect_equal(.ret$restarts, 0L)
  expect_lt(.ret$newtonDecrement, .control$outerTrustFterm)

  # An active bound: the unconstrained minimum (0, 0) is outside the box, so
  # every Newton step leaves it.  Rejecting those trials alone converges only
  # linearly onto the bound; the coordinates are held on it instead and the
  # bounded minimum (1, 1) is reached exactly, with the gradient still pointing
  # out (a KKT point, not a stationary one).
  .n <- 0L
  .ret <- .trustOuter(
    c(3, 2),
    fn = function(x) {
      .n <<- .n + 1L
      x[1]^2 + 4 * x[2]^2
    },
    gr = function(x) c(2 * x[1], 8 * x[2]),
    lower = c(1, 1),
    upper = c(Inf, Inf),
    control = .control
  )
  expect_equal(.ret$x, c(1, 1))
  expect_equal(.ret$convergence, 0L)
  expect_identical(.ret$activeBounds, 1:2)
  expect_equal(.ret$gradient, c(2, 8))
  expect_lt(.n, 15L)
  # only x2's bound is active here: the Newton direction is toward (2, 0), so
  # x1 walks 3 -> 2 and never reaches its own bound.  `activeChanges` counts the
  # hold, and nothing is released -- the release itself is covered below.
  .ret <- .trustOuter(
    c(3, 2),
    fn = function(x) (x[1] - 2)^2 + 4 * x[2]^2,
    gr = function(x) c(2 * (x[1] - 2), 8 * x[2]),
    lower = c(1, 1),
    upper = c(Inf, Inf),
    control = .control
  )
  expect_equal(.ret$x, c(2, 1), tolerance = 1e-6)
  expect_identical(.ret$activeBounds, 2L)
  expect_identical(.ret$activeChanges, 1L)
})

test_that("a bound the minimum does not sit on is held and then released", {
  # Rosenbrock from (-1.2, 1) with x2 >= 0: the first steps run down the valley
  # and leave the box through x2, which is held on 0.  With x1 free the gradient
  # at the converged point points back up, so the hold is released and the true
  # interior minimum (1, 1) is reached.  Asserting the mechanism, not just the
  # answer: `activeChanges` is 2 (one hold, one release) and nothing stays held.
  # Without the release this stalls on the bound with x2 == 0.
  .control <- list(
    fast = TRUE, sigdig = 3, maxOuterIterations = 200L,
    outerTrustHessian = "analytic", outerTrustRinit = 0.4,
    outerTrustRmax = 3.2, outerTrustRestarts = 2L,
    outerTrustFterm = 1e-11, outerTrustMterm = 1e-11,
    hessian = function(x, relStep) {
      matrix(c(1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
               -400 * x[1], 200), 2, 2)
    }
  )
  .ret <- .trustOuter(
    c(-1.2, 1),
    fn = function(x) 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2,
    gr = function(x) {
      c(-400 * x[1] * (x[2] - x[1]^2) - 2 * (1 - x[1]), 200 * (x[2] - x[1]^2))
    },
    lower = c(-Inf, 0),
    upper = c(Inf, Inf),
    control = .control
  )
  expect_equal(.ret$x, c(1, 1), tolerance = 1e-5)
  expect_equal(.ret$convergence, 0L)
  expect_identical(.ret$activeBounds, integer(0))
  expect_identical(.ret$activeChanges, 2L)
})

test_that("the outer gradient and Hessian are skipped on a trial that cannot be accepted", {
  # The saving is the whole point of the lazy path, and no assertion on the
  # answer can see it: a build that evaluated gr() and the Hessian at every
  # trial would return exactly the same optimum.  So count the callbacks --
  # a trial worse than the incumbent must cost one fn() and nothing else.
  .fn <- 0L
  .gr <- 0L
  .hess <- 0L
  .control <- list(
    fast = TRUE, sigdig = 3, maxOuterIterations = 100L,
    outerTrustHessian = "analytic", outerTrustRinit = 0.4,
    outerTrustRmax = 3.2, outerTrustRestarts = 0L,
    outerTrustFterm = 1e-11, outerTrustMterm = 1e-11,
    hessian = function(x, relStep) {
      .hess <<- .hess + 1L
      matrix(c(1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
               -400 * x[1], 200), 2, 2)
    }
  )
  .ret <- .trustOuter(
    c(-1.2, 1),
    fn = function(x) {
      .fn <<- .fn + 1L
      100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
    },
    gr = function(x) {
      .gr <<- .gr + 1L
      c(-400 * x[1] * (x[2] - x[1]^2) - 2 * (1 - x[1]), 200 * (x[2] - x[1]^2))
    },
    lower = c(-Inf, -Inf),
    upper = c(Inf, Inf),
    control = .control
  )
  expect_equal(.ret$x, c(1, 1), tolerance = 1e-5)
  # Rosenbrock from this start is rejection-heavy, so the skip has to bite
  expect_lt(.gr, .fn)
  expect_identical(.hess, .gr)
  expect_identical(.hess, .ret$hessianEvaluations)
})

test_that("a hold cut short by the iteration budget still returns a point in the box", {
  # The hold is signalled by a condition raised from objfun, so trust returns
  # through its error path, where `argument` is the trial that left the box and
  # not a point it ever accepted.  Taking that as the answer put a parameter
  # outside its bound; the incumbent is given back instead.  Every budget is
  # swept because which one stops mid-hold is not obvious from the outside.
  .control <- list(
    fast = TRUE, sigdig = 3, outerTrustRestarts = 0L,
    outerTrustFterm = 1e-11, outerTrustMterm = 1e-11,
    hessian = function(x, relStep) diag(c(2, 8))
  )
  .lower <- c(1, 1)
  .warn <- character()
  for (.it in seq_len(20)) {
    .control$maxOuterIterations <- .it
    .ret <- withCallingHandlers(
      .trustOuter(
        c(3, 2),
        fn = function(x) x[1]^2 + 4 * x[2]^2,
        gr = function(x) c(2 * x[1], 8 * x[2]),
        lower = .lower,
        upper = c(Inf, Inf),
        control = .control
      ),
      warning = function(w) {
        .warn <<- c(.warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_true(all(.ret$x >= .lower),
                info = paste("maxOuterIterations =", .it))
  }
  # RcppTrust reports the condition as "error in first/last call to objfun";
  # that is this driver's own control flow and must not reach the fit's runInfo
  expect_false(any(grepl("call to objfun", .warn, fixed = TRUE)))
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
    foceiControl(print = 0L, calcTables = FALSE, covMethod = "", outerOpt = "trust", ...)
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
  fitN <- .nlmixr(
    model,
    d,
    "focei",
    foceiControl(print = 0L, calcTables = FALSE, covMethod = "", outerOpt = "nlminb", fast = TRUE)
  )
  expect_equal(fitA$objf, fitN$objf, tolerance = 1e-3)
  expect_equal(unname(fixef(fitA)), unname(fixef(fitN)), tolerance = 1e-2)
})

test_that("outerOpt='trust' holds a coordinate on an active bound and converges", {
  skip_on_cran()
  # An upper bound below the unconstrained optimum: every Newton step from
  # inside the box points out through it.  Rejecting those trials only shrinks
  # the region, which converges linearly onto the bound and can stop there
  # with a large gradient.  The driver holds the coordinate on the bound
  # instead and finishes the others; the bounded optimum is the one L-BFGS-B
  # finds.
  model <- function() {
    ini({ tka <- c(-Inf, 0.1, 0.3); tcl <- 1; tv <- 3.45
          eta.cl ~ 0.3; add.sd <- 0.7 })
    model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
            d/dt(depot) <- -ka * depot
            d/dt(center) <- ka * depot - cl / v * center
            cp <- center / v
            cp ~ add(add.sd) })
  }
  d <- nlmixr2data::theo_sd
  ctl <- function(...) {
    foceiControl(print = 0L, calcTables = FALSE, covMethod = "", fast = TRUE, ...)
  }
  fitT <- .nlmixr(model, d, "focei", ctl(outerOpt = "trust"))
  fitL <- .nlmixr(model, d, "focei", ctl(outerOpt = "lbfgsb3c"))
  expect_equal(unname(fixef(fitT)["tka"]), 0.3, tolerance = 1e-4)
  expect_identical(fitT$env$optReturn$activeBounds, 1L)
  expect_true(fitT$env$optReturn$converged)
  expect_equal(fitT$objf, fitL$objf, tolerance = 1e-3)
  expect_equal(unname(fixef(fitT)), unname(fixef(fitL)), tolerance = 1e-2)
})
