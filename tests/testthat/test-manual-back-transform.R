test_that("manual back-transform", {

  t100 <- function(x) {
    x * 100
  }

  one.cmt <- function() {
    ini({
      tka <- fix(0.45); backTransform("t100")
      tcl <- log(c(0, 2.7, 100)); backTransform("t100")
      tv <- 3.45; backTransform("none")
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

  skip_if_not(rxode2::.linCmtSensB())

  fit <- nlmixr(one.cmt, theo_sd, est="saem", control=saemControl(print=0, nBurn = 1, nEm = 1))

  expect_equal(setNames(fit$parFixedDf["tka", "Estimate"] * 100, NULL),
               setNames(fit$parFixedDf["tka", "Back-transformed"], NULL))

  expect_equal(fit$parFixed["tka", "Back-transformed(95%CI)"],
               sprintf("%3g", fit$parFixedDf["tka", "Estimate"] * 100))

  expect_equal(NA_real_,
               setNames(fit$parFixedDf["tka", "CI Lower"], NULL))
  expect_equal(NA_real_,
               setNames(fit$parFixedDf["tka", "CI Upper"], NULL))

  expect_equal(setNames(fit$parFixedDf["tcl", "Estimate"] * 100, NULL),
               setNames(fit$parFixedDf["tcl", "Back-transformed"], NULL))


  qn <- qnorm(1.0 - (1 - 0.95) / 2)

  expect_equal(setNames(t100(fit$parFixedDf["tcl", "Estimate"] + qn * fit$parFixedDf["tcl", "SE"]), NULL),
               setNames(fit$parFixedDf["tcl", "CI Upper"], NULL))

  expect_equal(setNames(t100(fit$parFixedDf["tcl", "Estimate"] - qn * fit$parFixedDf["tcl", "SE"]), NULL),
               setNames(fit$parFixedDf["tcl", "CI Lower"], NULL))

  expect_equal(fit$parFixed["tcl", "Back-transformed(95%CI)"],
               sprintf("%3g (%3g, %3g)",
                       t100(fit$parFixedDf["tcl", "Estimate"]),
                       setNames(t100(fit$parFixedDf["tcl", "Estimate"] - qn * fit$parFixedDf["tcl", "SE"]), NULL),
                       setNames(t100(fit$parFixedDf["tcl", "Estimate"] + qn * fit$parFixedDf["tcl", "SE"]), NULL)))


  expect_equal(setNames(exp(fit$parFixedDf["tv", "Estimate"] + qn * fit$parFixedDf["tv", "SE"]), NULL),
               setNames(fit$parFixedDf["tv", "CI Upper"], NULL))

  expect_equal(setNames(exp(fit$parFixedDf["tv", "Estimate"] - qn * fit$parFixedDf["tv", "SE"]), NULL),
               setNames(fit$parFixedDf["tv", "CI Lower"], NULL))

  qn <- qnorm(1.0 - (1 - 0.80) / 2)

  fit <- nlmixr(one.cmt, theo_sd, est="saem", control=saemControl(print=0, nBurn = 1, nEm = 1, ci=0.8))

  expect_equal(setNames(fit$parFixedDf["tka", "Estimate"] * 100, NULL),
               setNames(fit$parFixedDf["tka", "Back-transformed"], NULL))

  expect_equal(fit$parFixed["tka", "Back-transformed(80%CI)"],
               sprintf("%3g", fit$parFixedDf["tka", "Estimate"] * 100))

  expect_equal(NA_real_,
               setNames(fit$parFixedDf["tka", "CI Lower"], NULL))
  expect_equal(NA_real_,
               setNames(fit$parFixedDf["tka", "CI Upper"], NULL))

  expect_equal(setNames(fit$parFixedDf["tcl", "Estimate"] * 100, NULL),
               setNames(fit$parFixedDf["tcl", "Back-transformed"], NULL))


  expect_equal(setNames(t100(fit$parFixedDf["tcl", "Estimate"] + qn * fit$parFixedDf["tcl", "SE"]), NULL),
               setNames(fit$parFixedDf["tcl", "CI Upper"], NULL))

  expect_equal(setNames(t100(fit$parFixedDf["tcl", "Estimate"] - qn * fit$parFixedDf["tcl", "SE"]), NULL),
               setNames(fit$parFixedDf["tcl", "CI Lower"], NULL))

  expect_equal(fit$parFixed["tcl", "Back-transformed(80%CI)"],
               sprintf("%3g (%3g, %3g)",
                       t100(fit$parFixedDf["tcl", "Estimate"]),
                       setNames(t100(fit$parFixedDf["tcl", "Estimate"] - qn * fit$parFixedDf["tcl", "SE"]), NULL),
                       setNames(t100(fit$parFixedDf["tcl", "Estimate"] + qn * fit$parFixedDf["tcl", "SE"]), NULL)))


  expect_equal(setNames(exp(fit$parFixedDf["tv", "Estimate"] + qn * fit$parFixedDf["tv", "SE"]), NULL),
               setNames(fit$parFixedDf["tv", "CI Upper"], NULL))

  expect_equal(setNames(exp(fit$parFixedDf["tv", "Estimate"] - qn * fit$parFixedDf["tv", "SE"]), NULL),
               setNames(fit$parFixedDf["tv", "CI Lower"], NULL))

  expect_error(nlmixr(one.cmt, theo_sd, est="focei",
                      control=foceiControl(print=0, maxOuterIterations = 0, maxInnerIterations = 0,
                                           covMethod = "")), NA)

})
