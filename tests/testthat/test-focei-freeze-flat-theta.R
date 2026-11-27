test_that("freezeFlatTheta fixes a theta that cannot move the objective (#1382)", {
  skip_on_cran()
  dead <- function() {
    ini({
      tka <- 0.5
      tcl <- -3
      tv <- 3
      tdead <- 1
      eta.cl ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      junk <- tdead * 2
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  fit <- suppressMessages(nlmixr2(
    dead,
    nlmixr2data::theo_sd,
    "focei",
    control = foceiControl(print = 0, freezeFlatTheta = TRUE, calcTables = FALSE)
  ))
  expect_true(any(grepl("'tdead' has a flat gradient", fit$runInfo, fixed = TRUE)))
  expect_equal(unname(fit$theta["tdead"]), 1)
  expect_true(fit$iniDf$fix[fit$iniDf$name == "tdead"])
  expect_false("tdead" %in% colnames(fit$cov))
  expect_true(abs(fit$theta["tcl"] - (-3)) > 0.5)
  expect_true(abs(fit$theta["tka"] - 0.5) > 0.01)

  withr::with_options(list(nlmixr2.retryFocei = FALSE), {
    fit2 <- suppressMessages(nlmixr2(
      dead,
      nlmixr2data::theo_sd,
      "focei",
      control = foceiControl(print = 0, freezeFlatTheta = TRUE, calcTables = FALSE)
    ))
  })
  expect_equal(unname(fit2$theta["tdead"]), 1)
  expect_true(fit2$iniDf$fix[fit2$iniDf$name == "tdead"])
})

test_that("freezeFlatTheta leaves a model whose thetas are all used alone", {
  skip_on_cran()
  one <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  fit <- suppressMessages(nlmixr2(
    one,
    nlmixr2data::theo_sd,
    "focei",
    control = foceiControl(print = 0, freezeFlatTheta = TRUE, calcTables = FALSE)
  ))
  expect_false(any(grepl("flat gradient", fit$runInfo, fixed = TRUE)))
  expect_false(any(fit$iniDf$fix))
})

test_that("freezeFlatTheta is validated", {
  expect_false(foceiControl()$freezeFlatTheta)
  expect_error(foceiControl(freezeFlatTheta = NA))
})
