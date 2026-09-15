test_that("detHessian='conditional' validates its control", {
  expect_error(foceiControl(detHessian = "conditional"), "requires fast")
  expect_error(foceiControl(fast = TRUE, interaction = FALSE, detHessian = "conditional"),
               "requires fast")
  expect_no_error(foceiControl(fast = TRUE, detHessian = "conditional"))
  expect_identical(foceiControl(fast = TRUE, detHessian = "conditional")$detHessian, "conditional")
  expect_identical(foceiControl()$detHessian, "focei")
})

test_that("the full conditional determinant matches its closed form and its gradient", {
  skip_on_cran()
  model <- function() {
    ini({ level <- 0.2; error <- fix(0.2); etaLevel ~ 0.2 })
    model({ prediction <- exp(level+etaLevel); prediction ~ prop(error) })
  }
  data <- data.frame(ID = rep(1:3, each = 3), TIME = rep(1:3, 3),
    DV = c(1.5,1.7,1.6,0.8,1,1.1,2,2.1,1.9), AMT = 0, EVID = 0)
  .ctl <- function(det) foceiControl(fast = TRUE, detHessian = det,
    maxOuterIterations = 0L, maxInnerIterations = 1000L, epsilon = 1e-10,
    trustFterm = 1e-12, trustMterm = 1e-12, print = 0, covMethod = "",
    calcTables = FALSE, compress = FALSE)
  fitFocei <- .nlmixr(model, data, "focei", control = .ctl("focei"))
  fitFull <- .nlmixr(model, data, "focei", control = .ctl("conditional"))
  # the objective differs from FOCEI only through the log-determinant
  expect_false(isTRUE(all.equal(fitFocei$objf, fitFull$objf)))
  expect_equal(fitFocei$eta, fitFull$eta, tolerance = 1e-5)
  # -2LL = sum_i 2*obj_i(eta*) + log|d2 obj_i/d eta2| + log|Omega|, up to the 2*pi constant
  reference <- vapply(split(data, data$ID), function(rows) {
    objective <- function(eta) {
      prediction <- exp(0.2+eta)
      -sum(dnorm(rows$DV, prediction, 0.2*prediction, log = TRUE)) + 0.5*eta^2/0.2
    }
    eta <- optimize(objective, c(-2, 2), tol = 1e-12)$minimum
    2*objective(eta) + log(numDeriv::hessian(objective, eta)[1, 1]) + log(0.2)
  }, numeric(1))
  expect_equal(fitFull$objf, sum(reference) - nrow(data)*log(2*pi), tolerance = 1e-6)
  # the analytic gradient carries the determinant's third-order term
  expect_gt(fitFull$env$nAnalyticGradDirect, 0)
  gradient <- .foceiGradDirect(fitFull)
  expect_false(is.null(gradient))
  objAt <- function(level) {
    .nlmixr(rxode2::ini(fitFull$finalUi, level = level), data, "focei",
            control = .ctl("conditional"))$objf
  }
  h <- 1e-4
  fdLevel <- (objAt(0.2 + h) - objAt(0.2 - h)) / (2*h)
  expect_equal(unname(gradient[1]), fdLevel, tolerance = 1e-3)
  # ... and vanishes at the optimum of the full-determinant objective
  converged <- .nlmixr(model, data, "focei", control = foceiControl(fast = TRUE,
    detHessian = "conditional", print = 0, covMethod = "", calcTables = FALSE,
    compress = FALSE))
  expect_gt(converged$env$nAnalyticGradDirect, 0)
  expect_true(all(abs(.foceiGradDirect(converged)) < 0.05))
})

test_that("the full determinant sends the analytic covariance and outer Hessian to FD", {
  skip_on_cran()
  model <- function() {
    ini({ level <- 0.2; error <- fix(0.2); etaLevel ~ 0.2 })
    model({ prediction <- exp(level+etaLevel); prediction ~ prop(error) })
  }
  data <- data.frame(ID = rep(1:3, each = 3), TIME = rep(1:3, 3),
    DV = c(1.5,1.7,1.6,0.8,1,1.1,2,2.1,1.9), AMT = 0, EVID = 0)
  fit <- .nlmixr(model, data, "focei", control = foceiControl(fast = TRUE,
    detHessian = "conditional", outerOpt = "trust", outerTrustHessian = "analytic",
    print = 0, covMethod = "r", calcTables = FALSE, compress = FALSE))
  expect_false(grepl("analytic", fit$covMethod))
  expect_true(isTRUE(fit$env$optReturn$hessianFallback))
})
