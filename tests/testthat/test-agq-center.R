test_that("agqCenter control and objective labels", {
  expect_equal(agqControl()$agqCenter, "mode")
  expect_equal(agqControl(agqCenter = "start")$agqCenter, "start")
  expect_error(agqControl(agqCenter = "foo"))
  expect_equal(.agqRuleLabel(agqControl()), "")
  expect_equal(.agqRuleLabel(agqControl(agqRule = "stieltjes")), "AGQ-Stieltjes")
  expect_equal(.agqRuleLabel(agqControl(agqRule = "stieltjes", agqCenter = "prior")), "GQ-Stieltjes")
  expect_equal(.agqRuleLabel(agqControl(nAGQ = 3, agqCenter = "start")), "GQ3")
  expect_equal(.agqRuleLabel(agqControl(nAGQ = 3, agqRule = "legendre")), "AGQ-GL3")
})

test_that("grids not centred on the EBE", {
  skip_on_cran()
  p <- pump
  p$dv <- p$y
  p$time <- 0
  f <- function() {
    ini({
      t1 <- 1
      t2 <- 1
      t3 <- 1
      t4 <- 1
      eta1 ~ 1
    })
    model({
      if (group == 1) {
        lp <- t1 + t2 * logtstd
      } else {
        lp <- t3 + t4 * logtstd
      }
      lp <- lp + eta1
      lam <- exp(lp)
      y ~ dpois(lam)
    })
  }
  .fit <- function(..., maxOuterIterations = 0) {
    suppressMessages(suppressWarnings(nlmixr2(
      f,
      p,
      est = "agq",
      control = agqControl(maxOuterIterations = maxOuterIterations, print = 0, ...)
    )))
  }
  .ref <- .fit(nAGQ = 9, covMethod = "")$objective
  # placed at the start, the grid is the adaptive grid of that evaluation
  .st <- .fit(agqRule = "stieltjes", covMethod = "")$objective
  expect_equal(.fit(agqRule = "stieltjes", agqCenter = "start", covMethod = "")$objective, .st)
  expect_equal(.fit(nAGQ = 9, agqCenter = "start", covMethod = "")$objective, .ref)
  # a prior grid converges to the same integral, slowly
  .pr <- .fit(nAGQ = 40, agqCenter = "prior", covMethod = "")
  expect_lt(abs(.pr$objective - .ref), 0.05)
  expect_equal(rownames(.pr$objDf), "GQ40")
  # an estimation run on a grid frozen at good estimates; no analytic gradient
  .fm <- .fit(agqRule = "stieltjes", maxOuterIterations = 1000)
  .fs <- suppressMessages(suppressWarnings(nlmixr2(
    .fm,
    p,
    est = "agq",
    control = agqControl(agqRule = "stieltjes", agqCenter = "start", print = 0)
  )))
  expect_equal(unname(fixef(.fs)[1:4]), unname(fixef(.fm)[1:4]), tolerance = 0.05)
  expect_true(is.finite(.fs$objective))
  expect_false(is.null(.fs$cov))
})
