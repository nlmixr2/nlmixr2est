test_that('outer-model CSE preserves covariate names in a fast fit', {
  model <- function() {
    ini({ level <- .211123; error <- .2; eta ~ .3 })
    model({ d/dt(dummy) <- 0; dummy(0) <- 1
            prediction <- exp(level + eta) * dummy + my_rx_expr_testcov
            prediction ~ add(error) })
  }
  data <- data.frame(ID = 1, TIME = 1:3, DV = c(1.6, 1.7, 1.8), my_rx_expr_testcov = c(.1, .2, .3), AMT = 0, EVID = 0)
  fit <- function(optimize) {
    suppressWarnings(suppressMessages(nlmixr2(
      model,
      data,
      'focei',
      foceiControl(
        fast = TRUE,
        optExpression = optimize,
        maxOuterIterations = 0,
        print = 0,
        covMethod = '',
        calcTables = FALSE
      )
    )))
  }
  reference <- fit(FALSE)
  optimized <- fit(TRUE)
  # each model goes through rxOptExpr() once; rxode2 warns on a second pass
  expect_false(any(grepl("optimized more than once", optimized$runInfo, fixed = TRUE)))
  expect_equal(as.numeric(optimized$objf), as.numeric(reference$objf), tolerance = 1e-8)
  expect_equal(as.matrix(optimized$eta[-1]), as.matrix(reference$eta[-1]), tolerance = 1e-8)
  expect_gt(optimized$env$nAnalyticGradDirect, 0)
  expect_equal(optimized$env$nFDGradFast, 0)
})
