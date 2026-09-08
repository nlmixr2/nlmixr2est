test_that("conditional curvature preserves the FOCEI marginal objective", {
  skip_on_cran()
  model <- function() {
    ini({ level <- 0.2; error <- fix(0.2); etaLevel ~ 0.2 })
    model({ prediction <- exp(level+etaLevel); prediction ~ prop(error) })
  }
  data <- data.frame(ID = rep(1:3, each = 3), TIME = rep(1:3, 3),
    DV = c(1.5,1.7,1.6,0.8,1,1.1,2,2.1,1.9), AMT = 0, EVID = 0)
  for (censoring in c("none","m2","m3","m4")) for (inner in c("trust", "n1qn1")) {
    fitData <- data
    if (censoring %in% c("m3","m4")) {
      fitData$CENS <- as.integer(fitData$DV <= 1)
      fitData$DV[fitData$CENS == 1] <- 1
    }
    if (censoring == "m2") fitData$CENS <- 0L
    if (censoring %in% c("m2","m4")) fitData$LIMIT <- 0
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
    reference <- vapply(split(fitData,fitData$ID), function(rows) {
      objective <- function(eta) {
        prediction <- exp(0.2+eta)
        sd <- 0.2*prediction
        ll <- dnorm(rows$DV,prediction,sd,log = TRUE)
        if (censoring %in% c("m2","m4")) ll <- ll-pnorm(0,prediction,sd,lower.tail = FALSE,log.p = TRUE)
        if (censoring %in% c("m3","m4")) {
          probability <- pnorm(rows$DV,prediction,sd)
          if (censoring == "m4") probability <- (probability-pnorm(0,prediction,sd))/pnorm(0,prediction,sd,lower.tail = FALSE)
          ll[rows$CENS == 1] <- log(probability[rows$CENS == 1])
        }
        -sum(ll)+0.5*eta^2/0.2
      }
      optimize(objective,c(-2,2),tol = 1e-10)$minimum
    }, numeric(1))
    expect_equal(as.numeric(fits[[2]]$eta[[2]]),unname(reference),tolerance = 1e-5)
  }
})

test_that("conditional curvature requires a compatible optimizer", {
  expect_error(foceiControl(innerHessian = "conditional"), "requires fast")
  expect_error(foceiControl(fast = TRUE, innerOpt = "BFGS",
    innerHessian = "conditional"), "requires fast")
  expect_error(foceiControl(fast = TRUE, interaction = FALSE,
    innerHessian = "conditional"), "requires fast")
})
