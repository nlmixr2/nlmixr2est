test_that("pooled outer curvature matches derivatives of the reported objective", {
  skip_on_cran()
  model <- function() {
    ini({ logCl <- 1; logV <- 3; error <- 0.4
          etaCl + etaV ~ c(0.2, 0.03, 0.15) })
    model({ cl <- exp(logCl+etaCl); v <- exp(logV+etaV)
            d/dt(central) <- -cl/v*central
            prediction <- central/v; prediction ~ add(error) })
  }
  data <- data.frame(ID = rep(1:2, each = 6), TIME = rep(c(0, 0.5, 1, 2, 4, 8), 2),
                     DV = c(NA, 4.7, 4.1, 3.5, 2.7, 1.7, NA, 4.4, 4, 3.1, 2, 0.9),
                     AMT = rep(c(100, rep(0, 5)), 2), EVID = rep(c(1, rep(0, 5)), 2), CMT = 1)
  for (inner in c("n1qn1", "trust")) for (objectiveScale in c(0, 100)) {
    observed <- NULL
    optimizer <- function(par, fn, gr, lower, upper, control) {
      fn(par); gradient <- gr(par); before <- fn(par)
      baseline <- fn(par)
      hessian <- control$hessian(par)
      after <- fn(par); afterGradient <- gr(par)
      halfStep <- control$hessian(par, relStep = 5e-4)
      reference <- matrix(0, length(par), length(par))
      for (j in seq_along(par)) {
        step <- 1e-3*max(1, abs(par[j]))
        plus <- minus <- par; plus[j] <- plus[j]+step; minus[j] <- minus[j]-step
        fn(plus); gp <- gr(plus)
        fn(minus); gm <- gr(minus)
        reference[, j] <- (gp-gm)/(2*step)
      }
      direction <- seq_along(par)/sqrt(sum(seq_along(par)^2))
      step <- 1e-2
      f0 <- fn(par)
      curvatureFd <- (fn(par+step*direction)-2*f0+fn(par-step*direction))/step^2
      observed <<- list(hessian = hessian, reference = (reference+t(reference))/2,
                        halfStep = halfStep, drift = abs(after-baseline),
                        normalDrift = abs(baseline-before), value = baseline,
                        gradientDrift = max(abs(afterGradient-gradient)),
                        curvature = drop(crossprod(direction, hessian %*% direction)),
                        curvatureFd = curvatureFd)
      fn(par)
      list(x = par, convergence = 0L, message = "Outer curvature check")
    }
    .nlmixr(model, data, "foceif", control = foceiControl(
      outerOpt = optimizer, innerOpt = inner, print = 0, covMethod = "", calcTables = FALSE,
      boundedTransform = FALSE, scaleObjective = objectiveScale, maxInnerIterations = 1000L,
      epsilon = 1e-10, trustFterm = 1e-12, trustMterm = 1e-12,
      rxControl = rxode2::rxControl(atol = 1e-11, rtol = 1e-11)))
    expect_type(observed, "list")
    expect_equal(observed$hessian, t(observed$hessian), tolerance = 1e-12)
    expect_lt(norm(observed$hessian-observed$reference, "F")/norm(observed$reference, "F"), 2e-3)
    expect_equal(observed$hessian, observed$halfStep, tolerance = 1e-5)
    expect_lte(observed$drift, max(1e-7*(1+abs(observed$value)), 2*observed$normalDrift))
    expect_lt(observed$gradientDrift, 1e-4)
    expect_lt(abs(observed$curvature-observed$curvatureFd)/max(1, abs(observed$curvature)), 5e-3)
  }
})
