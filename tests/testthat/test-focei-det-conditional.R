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

test_that("the full determinant's sigma and omega gradients match central differences", {
  skip_on_cran()
  model <- function() {
    ini({ level <- 0.2; error <- 0.2; etaLevel ~ 0.2 })
    model({ prediction <- exp(level+etaLevel); prediction ~ prop(error) })
  }
  data <- data.frame(ID = rep(1:3, each = 3), TIME = rep(1:3, 3),
    DV = c(1.5,1.7,1.6,0.8,1,1.1,2,2.1,1.9), AMT = 0, EVID = 0)
  .ctl <- foceiControl(fast = TRUE, detHessian = "conditional", diagXform = "sqrt",
    maxOuterIterations = 0L, maxInnerIterations = 1000L, epsilon = 1e-10,
    trustFterm = 1e-12, trustMterm = 1e-12, print = 0, covMethod = "",
    calcTables = FALSE, compress = FALSE)
  fit <- .nlmixr(model, data, "focei", control = .ctl)
  expect_gt(fit$env$nAnalyticGradDirect, 0)
  gradient <- .foceiGradDirect(fit)
  objAt <- function(ui) .nlmixr(ui, data, "focei", control = .ctl)$objf
  h <- 1e-4
  central <- function(up, down) (objAt(up) - objAt(down)) / (2*h)
  ui <- fit$finalUi
  expect_equal(gradient[["level"]],
               central(rxode2::ini(ui, level = 0.2 + h), rxode2::ini(ui, level = 0.2 - h)),
               tolerance = 1e-4)
  expect_equal(gradient[["error"]],
               central(rxode2::ini(ui, error = 0.2 + h), rxode2::ini(ui, error = 0.2 - h)),
               tolerance = 1e-4)
  # om.chol is var^(-1/4) under diagXform="sqrt", so d var/d om.chol = -4 var^(5/4)
  dVar <- central(eval(bquote(rxode2::ini(ui, etaLevel ~ .(0.2 + h)))),
                  eval(bquote(rxode2::ini(ui, etaLevel ~ .(0.2 - h)))))
  expect_equal(gradient[["om.chol.1"]], dVar * -4 * 0.2^1.25, tolerance = 1e-4)
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

test_that("the full determinant's gradient stays consistent with censored observations", {
  skip_on_cran()
  model <- function() {
    ini({ level <- 0.2; error <- fix(0.2); etaLevel ~ 0.2 })
    model({ prediction <- exp(level+etaLevel); prediction ~ prop(error) })
  }
  data <- data.frame(ID = rep(1:3, each = 3), TIME = rep(1:3, 3),
    DV = c(1.5,1.7,1.6,0.8,1,1.1,2,2.1,1.9), AMT = 0, EVID = 0,
    CENS = c(0,0,0,1,0,1,0,0,0))
  .ctl <- foceiControl(fast = TRUE, detHessian = "conditional",
    maxOuterIterations = 0L, maxInnerIterations = 1000L, epsilon = 1e-10,
    trustFterm = 1e-12, trustMterm = 1e-12, print = 0, covMethod = "",
    calcTables = FALSE, compress = FALSE)
  fit <- .nlmixr(model, data, "focei", control = .ctl)
  expect_gt(fit$env$nAnalyticGradDirect, 0)
  gradient <- .foceiGradDirect(fit)
  objAt <- function(level) {
    .nlmixr(rxode2::ini(fit$finalUi, level = level), data, "focei", control = .ctl)$objf
  }
  h <- 1e-4
  expect_equal(unname(gradient[1]), (objAt(0.2 + h) - objAt(0.2 - h)) / (2*h), tolerance = 1e-4)
})
