test_that("fast gradients use the reported objective scale", {
  skip_on_cran()
  model <- function() {
    ini({ level <- 0.2; error <- fix(0.5); etaLevel ~ 0.2 })
    model({ prediction <- exp(level+etaLevel); prediction ~ add(error) })
  }
  data <- data.frame(ID = c(1, 1, 2, 2), TIME = c(1, 2, 1, 2),
                     DV = c(1.5, 1.7, 0.8, 1), AMT = 0, EVID = 0)
  observed <- NULL
  optimizer <- function(par, fn, gr, lower, upper, control) {
    fn(par); gradient <- gr(par)
    reference <- vapply(seq_along(par), function(j) {
      shift <- rep(0, length(par)); shift[j] <- 1e-4
      (fn(par+shift)-fn(par-shift))/2e-4
    }, numeric(1))
    observed <<- list(gradient = gradient, reference = reference)
    fn(par)
    list(x = par, convergence = 0L, message = "Gradient scale check")
  }
  .nlmixr(model, data, "foceif", control = foceiControl(
    outerOpt = optimizer, innerOpt = "n1qn1", scaleObjective = 100,
    epsilon = 1e-10, maxInnerIterations = 1000L, covMethod = "", calcTables = FALSE, print = 0))
  expect_equal(observed$gradient, observed$reference, tolerance = 1e-4)
})
