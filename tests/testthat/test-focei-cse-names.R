nmTest({
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
})

nmTest({
  test_that("each fast FOCEi model is sum-prod'd and CSE-optimized only once", {
    seen <- new.env()
    seen$sp <- list(input = character(0), output = character(0))
    seen$oe <- list(input = character(0), output = character(0))
    track <- function(key, fun) {
      force(fun)
      function(x, ...) {
        txt <- paste(x, collapse = "\n")
        rec <- seen[[key]]
        # no model text is processed twice, and no processed output is fed back
        expect_false(txt %in% rec$input)
        expect_false(txt %in% rec$output)
        ret <- fun(x, ...)
        rec$input <- c(rec$input, txt)
        rec$output <- c(rec$output, paste(ret, collapse = "\n"))
        seen[[key]] <- rec
        ret
      }
    }
    local_mocked_bindings(
      rxSumProdModel = track("sp", rxode2::rxSumProdModel),
      rxOptExpr = track("oe", rxode2::rxOptExpr),
      .package = "rxode2"
    )
    model <- function() {
      ini({ tka <- log(1); tkel <- log(.2); error <- .2; eta.ka ~ .1 })
      model({ ka <- exp(tka + eta.ka); kel <- exp(tkel)
              d/dt(depot) <- -ka * depot
              d/dt(central) <- ka * depot - kel * central
              cp <- central / 10; cp ~ prop(error) })
    }
    data <- data.frame(
      ID = 1,
      TIME = c(0, 1, 2, 4),
      DV = c(NA, 6, 7, 5),
      AMT = c(100, 0, 0, 0),
      EVID = c(1, 0, 0, 0)
    )
    fit <- suppressWarnings(suppressMessages(nlmixr2(
      model,
      data,
      "focei",
      foceiControl(
        fast = TRUE,
        sumProd = TRUE,
        maxOuterIterations = 0,
        print = 0,
        covMethod = "",
        calcTables = FALSE
      )
    )))
    expect_gt(length(seen$sp$input), 0L)
    expect_gt(length(seen$oe$input), 0L)
    expect_gt(fit$env$nAnalyticGradDirect, 0)
  })
})
