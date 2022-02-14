test_that("testing saem without table can add focei objf", {

  one.cmt <- function() {
    ini({
      tka <- 0.45 ; label("Log Ka")
      tcl <- 1 ; label("Log Cl")
      tv <- 3.45 ; label("log V")
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  fit <-
    suppressMessages(
      nlmixr(one.cmt, theo_sd, est = "saem", control = list(calcTables = FALSE, print = 0))
    )

  expect_s3_class(fit, "nlmixr2FitCore")
  expect_false(inherits(fit, "data.frame"))
  expect_false(inherits(fit, "nlmixrFitData"))

  expect_error(suppressMessages(setOfv(fit, "focei")), NA)
  expect_error(suppressMessages(setOfv(fit, "foce")), NA)
  expect_error(suppressMessages(setOfv(fit, "fo")), NA)
  expect_error(suppressMessages(setOfv(fit, "gauss3_1.6")), NA)
  expect_error(suppressMessages(setOfv(fit, "laplace2")), NA)
  expect_error(suppressMessages(setOfv(fit, "laplace3")), NA)
})
