test_that("stieltjesControl defaults and conversions", {
  .c <- stieltjesControl()
  expect_s3_class(.c, "stieltjesControl")
  expect_equal(.c$agqRule, "stieltjes")
  expect_equal(.c$agqCenter, "start")
  expect_equal(.c$stage1, "focei")
  expect_error(stieltjesControl(nAGQ = 1), "nAGQ > 1")
  expect_error(stieltjesControl(stage1 = "saem"))
  # foceiControl-family controls convert
  .v <- suppressMessages(getValidNlmixrCtl.stieltjes(list(foceiControl(sigdig = 4))))
  expect_s3_class(.v, "stieltjesControl")
  expect_equal(.v$sigdig, 4)
  # the first stage runs its own method, without a covariance step
  .s1 <- .stieltjesStage1Control(stieltjesControl(stage1 = "foce", sigdig = 4))
  expect_s3_class(.s1, "foceControl")
  expect_equal(.s1$interaction, 0L)
  expect_equal(.s1$nAGQ, 0L)
  expect_equal(.s1$sigdig, 4)
  expect_equal(.s1$covMethod, 0L)
  expect_null(.s1$stage1)
  # deparse round trip
  .d <- rxode2::rxUiDeparse(stieltjesControl(stage1 = "laplace", agqGrid = c(4, 2, 0.5, 0.95)), "x")
  .x <- eval(.d)
  expect_equal(.x$stage1, "laplace")
  expect_equal(unname(.x$agqGrid), c(4, 2, 0.5, 0.95))
})

test_that("stage-1 estimates reach the internal bounded scale", {
  .ui <- new.env(parent = emptyenv())
  .ui$iniDf <- data.frame(
    name = c("ta", "rxBoundedTr.fr", "rxBoundedTr.sd", "eta.a"),
    neta1 = c(NA, NA, NA, 1),
    neta2 = c(NA, NA, NA, 1),
    est = c(0, 0, 0, 0.1)
  )
  .fit <- list(
    theta = c(ta = 2, fr = 0.25, sd = 0.5),
    omega = matrix(0.3, 1, 1)
  )
  .tr <- list(
    list(name = "fr", internalName = "rxBoundedTr.fr", type = "logit", lower = 0, upper = 1),
    list(name = "sd", internalName = "rxBoundedTr.sd", type = "lower_exp", lower = 0, upper = Inf)
  )
  .fit$omega <- matrix(0.3, 1, 1, dimnames = list("eta.a", "eta.a"))
  .stieltjesStage1Ini(.ui, .fit, .tr)
  expect_equal(.ui$iniDf$est, c(2, qlogis(0.25), log(0.5), 0.3))
  # an eta the stage-1 omega does not list keeps its value
  .ui$iniDf$est[4] <- 0.1
  .fit$omega <- matrix(0.3, 1, 1, dimnames = list("eta.z", "eta.z"))
  .stieltjesStage1Ini(.ui, .fit, .tr)
  expect_equal(.ui$iniDf$est[4], 0.1)
})

test_that("est='stieltjes' keeps the exact estimates of a model linear in eta", {
  skip_on_cran()
  set.seed(42)
  d <- do.call(
    rbind,
    lapply(1:20, function(i) {
      e <- rnorm(2, 0, c(0.7, 0.5))
      t <- c(1, 2, 4, 8)
      data.frame(id = i, time = t, dv = 1 + e[1] + (0.5 + e[2]) * t + rnorm(4, 0, 0.3))
    })
  )
  f <- function() {
    ini({
      ta <- 1
      tb <- 0.5
      eta.a ~ 0.5
      eta.b ~ 0.2
      add.sd <- 0.3
    })
    model({
      a <- ta + eta.a
      b <- tb + eta.b
      y <- a + b * time
      y ~ add(add.sd)
    })
  }
  .fit <- suppressMessages(suppressWarnings(nlmixr2(
    f,
    d,
    est = "stieltjes",
    control = stieltjesControl(print = 0)
  )))
  expect_equal(.fit$est, "stieltjes")
  expect_equal(rownames(.fit$objDf), "GQ-Stieltjes")
  expect_s3_class(getValidNlmixrCtl.stieltjes(list(.fit$control)), "stieltjesControl")
  .s1 <- .fit$stieltjesStage1
  expect_equal(.s1$est, "focei")
  # FOCEi is exact here, and the second stage's grid is exact where it starts
  expect_equal(.fit$objective, .s1$objective, tolerance = 1e-3)
  expect_equal(unname(.fit$theta[1:3]), unname(.s1$theta[1:3]), tolerance = 1e-2)
  # a fixed omega is carried through both stages
  .fx <- suppressMessages(suppressWarnings(nlmixr2(
    rxode2::ini(f, eta.b = fix(0.2)),
    d,
    est = "stieltjes",
    control = stieltjesControl(print = 0, covMethod = "")
  )))
  expect_equal(unname(.fx$omega["eta.b", "eta.b"]), 0.2)
  expect_true(is.finite(.fx$objective))
})
