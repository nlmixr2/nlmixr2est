# covMethod = "auto" (#1217): analytic or "r,s" chosen per model from its size

test_that("foceiControl() defaults to covMethod = \"auto\" and round-trips it", {
  .ctl <- foceiControl()
  expect_identical(.ctl$covType, "auto")
  expect_identical(.ctl$covMethod, 1L)
  expect_identical(.foceiControlCovMethodName(.ctl), "auto")
  expect_identical(.foceiControlCovMethodName(foceiControl(covMethod = "r,s")), "r,s")
  expect_identical(.foceiControlCovMethodName(foceiControl(covMethod = "analytic")), "analytic")
  .re <- do.call(foceiControl, .ctl)
  expect_identical(.re$covType, "auto")
  # a resolved control still names the request it came from
  .res <- .ctl
  .res$covType <- "analytic"
  .res$covMethod <- 2L
  .res$covAuto <- TRUE
  expect_identical(.foceiControlCovMethodName(.res), "auto")
  expect_true(do.call(foceiControl, .res)$covAuto)
})

test_that("the predicted analytic / r,s cost ratio rises with states and falls with parameters", {
  expect_gt(.covAutoRatio(60, 1, 4), .covAutoRatio(2, 1, 4))
  expect_lt(.covAutoRatio(5, 8, 9), .covAutoRatio(5, 1, 9))
  expect_lt(.covAutoRatio(5, 3, 20), .covAutoRatio(5, 3, 4))
  expect_lt(.covAutoRatio(2, 3, 5), 1)
})

.covAutoOde <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

test_that(".foceiCovAuto resolves the request per model", {
  .ui <- rxode2::rxode2(.covAutoOde)
  .r <- .foceiCovAuto(.ui, foceiControl())
  expect_identical(.r$covType, "analytic")
  expect_identical(.r$covMethod, 2L)
  expect_true(.r$covAuto)
  # an explicit request is left alone
  .r <- .foceiCovAuto(.ui, foceiControl(covMethod = "r,s"))
  expect_identical(.r$covType, "fd")
  expect_false(.r$covAuto)
  # FO, linCmt() and a model the cost rule sends to r,s take the finite-difference step
  .fo <- foceiControl()
  class(.fo) <- "foControl"
  expect_identical(.foceiCovAuto(.ui, .fo)$covType, "fd")
  .lin <- rxode2::rxode2(function() {
    ini({
      tcl <- 1
      tv <- 3.45
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  })
  expect_identical(.foceiCovAuto(.lin, foceiControl())$covType, "fd")
  local_mocked_bindings(.covAutoRatio = function(...) 2)
  .r <- .foceiCovAuto(.ui, foceiControl())
  expect_identical(.r$covType, "fd")
  expect_identical(.r$covMethod, 1L)
  expect_true(.r$covAuto)
})

test_that("covMethod = \"auto\" installs the analytic covariance, or r,s quietly when it declines", {
  skip_on_cran()
  .fit <- suppressMessages(nlmixr2(.covAutoOde, theo_sd, "focei", control = foceiControl(print = 0L)))
  expect_identical(.fit$covMethod, "analytic (full)")
  expect_identical(.fit$env$covAuto, "analytic")
  # cholSECov is out of analytic scope, found only at the covariance step
  .msg <- character(0)
  .fit2 <- withCallingHandlers(
    nlmixr2(.covAutoOde, theo_sd, "focei", control = foceiControl(print = 0L, cholSECov = TRUE)),
    message = function(m) {
      .msg <<- c(.msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_false(any(grepl("analytic", .msg)))
  expect_match(.fit2$covMethod, "^r,s")
  expect_match(.fit2$env$covAuto, "^r,s: analytic declined \\(cholSECov")
  expect_false(any(grepl("analytic", .fit2$runInfo)))
  expect_false(any(.fit2$env$covTried$method == "analytic"))
  .rs <- suppressMessages(nlmixr2(
    .covAutoOde,
    theo_sd,
    "focei",
    control = foceiControl(print = 0L, cholSECov = TRUE, covMethod = "r,s")
  ))
  expect_identical(.fit2$covMethod, .rs$covMethod)
  expect_equal(.fit2$cov, .rs$cov, tolerance = 1e-6)
  expect_null(.rs$env$covAuto)
  local_mocked_bindings(.covAutoRatio = function(...) 2)
  .fit3 <- suppressMessages(nlmixr2(.covAutoOde, theo_sd, "focei", control = foceiControl(print = 0L)))
  expect_match(.fit3$covMethod, "^r,s")
  expect_match(.fit3$env$covAuto, "^r,s: out of analytic scope or predicted")
})
