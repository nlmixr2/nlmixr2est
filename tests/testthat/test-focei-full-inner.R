test_that("the f* sugar methods force the full conditional inner Hessian", {
  for (est in c("flaplace", "mflaplace", "iflaplace", "fagq", "mfagq", "ifagq")) {
    expect_true(est %in% nlmixr2AllEst(), info = est)
  }
  # NULL control: each f* validator has to reach its own base constructor's
  # defaults (nAGQ, muModel) on top of forcing the conditional curvature
  .ctl <- function(est) getValidNlmixrCtl(structure(list(NULL), class = est))
  for (est in c("flaplace", "mflaplace", "iflaplace", "fagq", "mfagq", "ifagq")) {
    .c <- .ctl(est)
    expect_true(.c$fast, info = est)
    expect_identical(.c$innerHessian, "conditional", info = est)
  }
  expect_equal(.ctl("flaplace")$nAGQ, 1L)
  expect_gt(.ctl("fagq")$nAGQ, 1)
  expect_identical(.ctl("mflaplace")$muModel, "lin")
  expect_identical(.ctl("iflaplace")$muModel, "irls")
  expect_identical(.ctl("mfagq")$muModel, "lin")
  expect_identical(.ctl("ifagq")$muModel, "irls")
})

test_that("a full conditional Laplace/AGQ fit says so", {
  skip_on_cran()
  model <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
          eta.cl ~ 0.3; add.sd <- 0.7 })
    model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv)
            d/dt(depot) <- -ka * depot
            d/dt(center) <- ka * depot - cl / v * center
            cp <- center / v
            cp ~ add(add.sd) })
  }
  d <- nlmixr2data::theo_sd
  .ctl <- function(...) {
    list(print = 0L, calcTables = FALSE, covMethod = "",
         maxOuterIterations = 0L, ...)
  }
  fitL <- .nlmixr(model, d, "flaplace", do.call(laplaceControl, .ctl()))
  expect_identical(fitL$method, "Full Laplace")
  expect_identical(rownames(fitL$objDf)[1], "Full Laplace")
  # broom/setOfv match the row by its lowercased label
  expect_identical(fitL$env$ofvType, tolower(rownames(fitL$objDf)[1]))
  expect_gt(fitL$env$nConditionalInnerHessian, 0L)

  fitA <- .nlmixr(model, d, "fagq", do.call(agqControl, .ctl(nAGQ = 3)))
  expect_identical(fitA$method, "Full AGQ")
  expect_identical(rownames(fitA$objDf)[1], "Full AGQ3")
  expect_identical(fitA$env$ofvType, tolower(rownames(fitA$objDf)[1]))
  expect_gt(fitA$env$nConditionalInnerHessian, 0L)

  # the base methods are untouched, and the conditional curvature does not move
  # the marginal objective it is curvature FOR
  base <- .nlmixr(model, d, "laplace", do.call(laplaceControl, .ctl()))
  expect_identical(base$method, "AGQ")
  expect_identical(rownames(base$objDf)[1], "Laplace")
  expect_identical(base$env$nConditionalInnerHessian, 0L)
  expect_equal(fitL$objf, base$objf, tolerance = 1e-4)
})

test_that("the full conditional inner Hessian refuses a generalized likelihood", {
  skip_on_cran()
  model <- function() {
    ini({ tf <- 0.5; eta.f ~ 0.1 })
    model({ p <- expit(tf + eta.f); n1 <- 1; dv ~ dbinom(n1, p) })
  }
  d <- data.frame(ID = rep(1:4, each = 3), TIME = rep(1:3, 4),
                  dv = c(1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0), AMT = 0, EVID = 0)
  expect_error(.nlmixr(model, d, "flaplace",
                       laplaceControl(print = 0L, calcTables = FALSE,
                                      covMethod = "", maxOuterIterations = 0L)),
               "Gaussian endpoints")
})
