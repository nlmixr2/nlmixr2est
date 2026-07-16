test_that("zero etas are restored in the final ui even with nested nlmixr2 calls (#741)", {
  skip_on_cran()

  set.seed(42)
  d <- do.call(rbind, lapply(1:12, function(i) {
    b <- 2 + rnorm(1, 0, 0.3)
    data.frame(id = i, time = 1:5, DV = 1 + b * (1:5) + rnorm(5, 0, 0.2))
  }))

  mod <- function() {
    ini({
      a <- 1
      b <- 2
      bsva ~ 0
      bsvb ~ 0.1
      addSd <- 1
    })
    model({
      ai <- a + bsva
      bi <- b + bsvb
      ci <- ai + bi * time
      ci ~ add(addSd)
    })
  }

  # est="posthoc" adds the focei objective through a nested nlmixr2() call
  # (.setOfvFo), which used to reset the zero-eta restore info before the
  # original model was put back, dropping bsva from the fit ui
  fit <- suppressMessages(suppressWarnings(
    nlmixr2(mod, data = d, est = "posthoc",
            control = foceiControl(print = 0))))

  expect_true("bsva" %in% fit$ui$iniDf$name)
  expect_true(grepl("bsva", paste(deparse(body(fit$ui$fun)), collapse = " ")))

  # and the eta can be re-enabled by piping
  newUi <- suppressMessages(rxode2::ini(fit, bsva ~ 0.1))
  expect_equal(newUi$iniDf$est[newUi$iniDf$name == "bsva"], 0.1)
})
