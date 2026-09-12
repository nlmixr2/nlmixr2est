test_that("built-in nlminb automatically uses supported outer curvature", {
  skip_on_cran()
  model <- function() {
    ini({ level <- 0.2; error <- fix(0.2); etaLevel ~ 0.2 })
    model({ prediction <- exp(level+etaLevel); prediction ~ prop(error) })
  }
  data <- data.frame(ID = rep(1:6, each = 3), TIME = rep(1:3, 6),
    DV = c(1.5, 1.7, 1.6, 0.8, 1, 1.1, 1.1, 1.2, 1, 2, 2.1, 1.9,
           0.9, 0.8, 0.7, 1.5, 1.3, 1.4), AMT = 0, EVID = 0)
  gradientOptimizer <- function(par, fn, gr, lower, upper, control) {
    control$hessian <- NULL
    .nlminb(par, fn, gr, lower, upper, control)
  }
  for (family in c("FOCEI", "FOCE+")) for (inner in c("n1qn1", "trust")) {
    control <- foceiControl(fast = TRUE, outerOpt = "nlminb", innerOpt = inner,
      epsilon = 1e-10, trustFterm = 1e-12, trustMterm = 1e-12,
      interaction = family == "FOCEI", foce = if (family == "FOCE+") "foce+" else "nonmem",
      maxInnerIterations = 1000L, maxOuterIterations = 1000L,
      rel.tol = 1e-8, x.tol = 1e-8, print = 0, covMethod = "", calcTables = FALSE, compress = FALSE)
    fit <- .nlmixr(model, data, "focei", control = control)
    control$outerOpt <- -1L; control$outerOptFun <- gradientOptimizer
    reference <- .nlmixr(model, data, "focei", control = control)
    expect_gt(fit$env$optReturn$hessianEvaluations, 0L)
    expect_false(fit$env$optReturn$hessianFallback)
    expect_equal(fit$env$optReturn$convergence, 0L)
    expect_identical(reference$env$optReturn$hessianEvaluations, 0L)
    expect_lt(abs(fit$objf-reference$objf), 1e-4)
    expect_equal(fit$theta, reference$theta, tolerance = 1e-3)
  }
  fallbackData <- data
  fallbackData$CENS <- as.integer(fallbackData$DV < 1)
  fallbackData$DV[fallbackData$CENS == 1] <- 1
  fallback <- .nlmixr(model, fallbackData, "focei", control = foceiControl(
    fast = TRUE, censOption = "laplace", outerOpt = "nlminb", innerOpt = "n1qn1",
    print = 0, covMethod = "", calcTables = FALSE, compress = FALSE))
  expect_true(fallback$env$optReturn$hessianFallback)
  expect_identical(fallback$env$optReturn$hessianEvaluations, 1L)
  expect_true(is.finite(fallback$objf))
})
