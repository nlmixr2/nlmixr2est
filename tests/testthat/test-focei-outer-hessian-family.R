test_that("FOCE and AGQ curvature includes their objective-specific terms", {
  skip_on_cran()
  model <- function() {
    ini({ level <- 0.2; slope <- -0.1; error <- 0.2
          etaLevel + etaSlope ~ c(0.2, 0.02, 0.05) })
    model({ prediction <- exp(level+etaLevel+(slope+etaSlope)*TIME)
            prediction ~ prop(error) })
  }
  data <- data.frame(ID = rep(1:2, each = 3), TIME = rep(1:3, 2),
    DV = c(1.5, 1.4, 1.3, 0.8, 0.85, 0.9), AMT = 0, EVID = 0)
  for (family in c("FOCE", "AGQ3", "AGQ5")) for (inner in c("n1qn1", "trust")) {
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
        reference[, j] <- (gp-gm)/4e-4
      }
      direction <- seq_along(par)/sqrt(sum(seq_along(par)^2))
      step <- 2e-3
      f0 <- fn(par); fp <- fn(par+step*direction); fm <- fn(par-step*direction)
      observed <<- list(hessian = hessian, reference = (reference+t(reference))/2,
        drift = abs(after-value), slope = sum(gradient*direction), slopeFd = (fp-fm)/(2*step),
        curvature = drop(crossprod(direction, hessian %*% direction)),
        curvatureFd = (fp-2*f0+fm)/step^2)
      fn(par)
      list(x = par, convergence = 0L, message = "FOCE/AGQ curvature check")
    }
    control <- foceiControl(fast = TRUE, outerOpt = optimizer, innerOpt = inner,
      interaction = family != "FOCE", nAGQ = if (family == "FOCE") 0L else as.integer(sub("AGQ", "", family)),
      print = 0, covMethod = "", calcTables = FALSE, boundedTransform = FALSE,
      epsilon = 1e-10, maxInnerIterations = 1000L, trustFterm = 1e-12, trustMterm = 1e-12,
      scaleType = "mult", scaleTo = 0,
      rxControl = rxode2::rxControl(atol = 1e-11, rtol = 1e-11))
    .nlmixr(model, data, "focei", control = control)
    expect_type(observed, "list")
    expect_equal(observed$hessian, t(observed$hessian), tolerance = 1e-10)
    expect_lt(norm(observed$hessian-observed$reference, "F")/norm(observed$reference, "F"), 2e-3)
    expect_lt(observed$drift, 1e-7)
    expect_equal(observed$slope, observed$slopeFd, tolerance = 2e-3)
    expect_equal(observed$curvature, observed$curvatureFd, tolerance = 2e-3)
  }
})
