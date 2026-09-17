test_that("a failed 3rd-order probe loosens the tolerance before falling back", {
  skip_on_cran()
  model <- function() {
    ini({ logCl <- 1; logV <- 3; error <- 0.4; etaCl ~ 0.2 })
    model({ cl <- exp(logCl+etaCl); v <- exp(logV)
            d/dt(central) <- -cl/v*central
            prediction <- central/v; prediction ~ add(error) })
  }
  data <- data.frame(ID = rep(1:2, each = 6), TIME = rep(c(0, 0.5, 1, 2, 4, 8), 2),
                     DV = c(NA, 4.7, 4.1, 3.5, 2.7, 1.7, NA, 4.4, 4, 3.1, 2, 0.9),
                     AMT = rep(c(100, rep(0, 5)), 2), EVID = rep(c(1, rep(0, 5)), 2), CMT = 1)
  # NLMIXR2EST_HESS_PROBE_FAIL fails the first N probe ATTEMPTS, which is the only
  # way to reach the state deterministically: no model-level knob makes a probe fail
  # at 1e-12 and solve a rung looser.
  .run <- function(fail) {
    Sys.setenv(NLMIXR2EST_HESS_PROBE_FAIL = fail)
    on.exit(Sys.unsetenv("NLMIXR2EST_HESS_PROBE_FAIL"), add = TRUE)
    .hessian <- NULL
    .optimizer <- function(par, fn, gr, lower, upper, control) {
      fn(par)
      .hessian <<- tryCatch(control$hessian(par), error = function(e) e)
      list(x = par, convergence = 0L, message = "probe tolerance check")
    }
    .fit <- .nlmixr(model, data, "focei", control = foceiControl(
      fast = TRUE, outerOpt = .optimizer, print = 0, covMethod = "",
      calcTables = FALSE, maxInnerIterations = 1000L, epsilon = 1e-10))
    list(hessian = .hessian, relax = .fit$env$nHessTolRelax)
  }
  # nothing fails -> the probes solve at the tightened tolerance, nothing loosened
  .clean <- .run("")
  expect_true(is.matrix(.clean$hessian))
  expect_equal(as.integer(.clean$relax), 0L)
  # one forced failure -> the rung looser rescues it, so the Hessian still arrives
  .retried <- .run("1")
  expect_true(is.matrix(.retried$hessian))
  expect_gt(as.integer(.retried$relax), 0L)
  expect_equal(.retried$hessian, .clean$hessian, tolerance = 1e-4)
  # every rung fails -> the budget is spent and only then does it give up
  .spent <- .run("99")
  expect_false(is.matrix(.spent$hessian))
  expect_gt(as.integer(.spent$relax), 0L)
})
