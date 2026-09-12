test_that("M2/M3/M4 outer curvature follows the censored fit objective", {
  skip_on_cran()
  makeModel <- function(errorType) {
    if (errorType == "add") function() {
      ini({ level <- 0.2; slope <- -0.1; error <- 0.5
            etaLevel + etaSlope ~ c(0.2, 0.02, 0.05) })
      model({ prediction <- exp(level+etaLevel+(slope+etaSlope)*TIME)
              prediction ~ add(error) })
    } else function() {
      ini({ level <- 0.2; slope <- -0.1; error <- 0.3
            etaLevel + etaSlope ~ c(0.2, 0.02, 0.05) })
      model({ prediction <- exp(level+etaLevel+(slope+etaSlope)*TIME)
              prediction ~ prop(error) })
    }
  }
  base <- data.frame(ID = rep(1:2, each = 5), TIME = rep(1:5, 2),
    DV = c(1.8, 1.6, 1.4, 1.2, 0.9, 0.9, 0.8, 0.75, 0.7, 0.6), AMT = 0, EVID = 0)
  for (errorType in c("add", "prop")) for (kind in c("M2", "M3", "M4"))
    for (family in c("FOCEI", "FOCE", "FOCE+", "AGQ")) {
      data <- base
      data$CENS <- 0L
      if (kind != "M2") {
        data$CENS <- ifelse(data$DV < 0.8, 1L, ifelse(data$DV > 1.7, -1L, 0L))
        data$DV[data$CENS == 1] <- 0.8; data$DV[data$CENS == -1] <- 1.7
      }
      if (kind != "M3") data$LIMIT <- ifelse(data$CENS == -1, 3, 0)
      observed <- NULL
      optimizer <- function(par, fn, gr, lower, upper, control) {
        fn(par); gradient <- gr(par); value <- fn(par)
        hessian <- control$hessian(par)
        after <- fn(par)
        reference <- matrix(0, length(par), length(par))
        for (j in seq_along(par)) {
          delta <- rep(0, length(par)); delta[j] <- 2e-4
          fn(par+delta); gp <- gr(par+delta)
          fn(par-delta); gm <- gr(par-delta)
          reference[,j] <- (gp-gm)/4e-4
        }
        direction <- seq_along(par)/sqrt(sum(seq_along(par)^2))
        step <- 2e-3
        f0 <- fn(par); fp <- fn(par+step*direction); fm <- fn(par-step*direction)
        observed <<- list(hessian = hessian, reference = (reference+t(reference))/2,
          drift = abs(after-value), slope = sum(gradient*direction), slopeFd = (fp-fm)/(2*step),
          curvature = drop(crossprod(direction, hessian %*% direction)),
          curvatureFd = (fp-2*f0+fm)/step^2)
        fn(par)
        list(x = par, convergence = 0L, message = "Censored curvature check")
      }
      .nlmixr(makeModel(errorType), data, "focei", control = foceiControl(
        fast = TRUE, outerOpt = optimizer, innerOpt = "n1qn1", censOption = "gauss",
        interaction = !(family %in% c("FOCE", "FOCE+")), foce = if (family == "FOCE+") "foce+" else "nonmem",
        nAGQ = if (family == "AGQ") 3L else 0L, print = 0, covMethod = "", calcTables = FALSE,
        boundedTransform = FALSE, epsilon = 1e-10, maxInnerIterations = 1000L,
        scaleType = "mult", scaleTo = 0,
        rxControl = rxode2::rxControl(atol = 1e-11, rtol = 1e-11)))
      expect_type(observed, "list")
      expect_equal(observed$hessian, t(observed$hessian), tolerance = 1e-10)
      expect_lt(norm(observed$hessian-observed$reference, "F")/norm(observed$reference, "F"), 3e-3)
      expect_lt(observed$drift, 1e-7)
      expect_equal(observed$slope, observed$slopeFd, tolerance = 3e-3)
      expect_equal(observed$curvature, observed$curvatureFd, tolerance = 3e-3)
    }
})

test_that("censored laplace keeps a consistent gradient and unavailable analytical curvature", {
  skip_on_cran()
  model <- function() {
    ini({ level <- 0.2; error <- fix(0.5); etaLevel ~ 0.2 })
    model({ prediction <- exp(level+etaLevel); prediction ~ add(error) })
  }
  data <- data.frame(ID = rep(1:2, each = 3), TIME = rep(1:3, 2),
    DV = c(1.5, 1.7, 1.6, 1, 1, 1.1), CENS = c(0, 0, 0, 1, 1, 0), AMT = 0, EVID = 0)
  for (family in c("FOCEI", "FOCE", "FOCE+", "AGQ")) {
    observed <- NULL
    optimizer <- function(par, fn, gr, lower, upper, control) {
      fn(par); gradient <- gr(par)
      hessian <- tryCatch(control$hessian(par), error = conditionMessage)
      reference <- vapply(seq_along(par), function(j) {
        delta <- rep(0, length(par)); delta[j] <- 1e-4
        (fn(par+delta)-fn(par-delta))/2e-4
      }, numeric(1))
      observed <<- list(gradient = gradient, reference = reference, hessian = hessian)
      fn(par)
      list(x = par, convergence = 0L, message = "Censored laplace scope check")
    }
    .nlmixr(model, data, "focei", control = foceiControl(
      fast = TRUE, outerOpt = optimizer, innerOpt = "n1qn1", censOption = "laplace",
      interaction = !(family %in% c("FOCE", "FOCE+")), foce = if (family == "FOCE+") "foce+" else "nonmem",
      nAGQ = if (family == "AGQ") 3L else 0L, print = 0, covMethod = "", calcTables = FALSE,
      epsilon = 1e-10, maxInnerIterations = 1000L))
    expect_equal(observed$gradient, observed$reference, tolerance = 3e-3)
    expect_match(observed$hessian, "status -4")
  }
})
