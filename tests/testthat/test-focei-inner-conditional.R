test_that("conditional curvature preserves the FOCEI marginal objective", {
  skip_on_cran()
  model <- function() {
    ini({ level <- 0.2; error <- fix(0.2); etaLevel ~ 0.2 })
    model({ prediction <- exp(level+etaLevel); prediction ~ prop(error) })
  }
  data <- data.frame(ID = rep(1:3, each = 3), TIME = rep(1:3, 3),
    DV = c(1.5,1.7,1.6,0.8,1,1.1,2,2.1,1.9), AMT = 0, EVID = 0)
  for (censored in c(FALSE,TRUE)) for (inner in c("trust", "n1qn1")) {
    fitData <- data
    if (censored) {
      fitData$CENS <- as.integer(fitData$DV <= 1)
      fitData$DV[fitData$CENS == 1] <- 1
    }
    fits <- lapply(c("focei", "conditional"), function(curvature) {
      .nlmixr(model, fitData, "focei", control = foceiControl(fast = TRUE,
        innerOpt = inner, innerHessian = curvature, maxOuterIterations = 0L,
        maxInnerIterations = 1000L, epsilon = 1e-10,
        trustFterm = 1e-12, trustMterm = 1e-12,
        print = 0, covMethod = "", calcTables = FALSE, compress = FALSE))
    })
    expect_equal(fits[[1]]$objf, fits[[2]]$objf, tolerance = 1e-6)
    expect_equal(fits[[1]]$eta, fits[[2]]$eta, tolerance = 1e-5)
    expect_identical(fits[[1]]$env$nConditionalInnerHessian, 0L)
    expect_gt(fits[[2]]$env$nConditionalInnerHessian, 0L)
  }
})

test_that("conditional curvature requires a compatible optimizer", {
  expect_error(foceiControl(innerHessian = "conditional"), "requires fast")
  expect_error(foceiControl(fast = TRUE, innerOpt = "BFGS",
    innerHessian = "conditional"), "requires fast")
  expect_error(foceiControl(fast = TRUE, interaction = FALSE,
    innerHessian = "conditional"), "requires fast")
})
