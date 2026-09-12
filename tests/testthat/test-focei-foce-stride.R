test_that("FOCE reads all ETA columns with the allocated observation stride", {
  skip_on_cran()
  model <- function() {
    ini({ level <- 0.2; slope <- -0.1; error <- fix(0.2)
          etaLevel ~ 0.2; etaSlope ~ 0.05 })
    model({ prediction <- exp(level+etaLevel+(slope+etaSlope)*TIME)
            prediction ~ prop(error) })
  }
  data <- data.frame(ID = rep(1:2, each = 3), TIME = rep(1:3, 2),
    DV = c(1.5, 1.4, 1.3, 0.8, 0.85, 0.9), AMT = 0, EVID = 0)
  fit <- .nlmixr(model, data, "foce", control = foceiControl(
    interaction = FALSE, innerOpt = "n1qn1", epsilon = 1e-10,
    maxOuterIterations = 0L, maxInnerIterations = 1000L, print = 0,
    covMethod = "", calcTables = FALSE, compress = FALSE))
  omega <- as.matrix(fit$omega); oi <- solve(omega)
  value <- 0
  for (id in 1:2) {
    eta <- as.numeric(fit$eta[id, -1]); time <- 1:3
    f <- exp(0.2+eta[1]+(-0.1+eta[2])*time)
    variance <- (0.2*exp(0.2-0.1*time))^2
    a <- cbind(f, time*f)
    hessian <- oi+crossprod(a, a/variance)
    value <- value+sum(log(variance)+(data$DV[data$ID == id]-f)^2/variance)+
      drop(crossprod(eta, oi %*% eta))+as.numeric(determinant(omega)$modulus)+
      as.numeric(determinant(hessian)$modulus)
  }
  expect_equal(fit$objf, value, tolerance = 1e-8)
})
