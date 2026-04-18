nmTest({
  skip_on_cran()

  # Shared model with one bounded parameter (td1 in [0, 1]).
  .logitModel <- function() {
    ini({
      tka   <- 0.45
      tcl   <- 1.0
      tv    <- 3.45
      td1   <- c(0, 0.5, 1)
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v  ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v  <- exp(tv  + eta.v)
      d1 <- td1
      f(depot) <- d1
      linCmt() ~ add(add.sd)
    })
  }

  # Lower-bound-only model (tlag >= 0).
  .lowerModel <- function() {
    ini({
      tka    <- log(1.57)
      tcl    <- log(0.0615)
      tv     <- log(3.5)
      tlag   <- c(0, 0.5)
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v  ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka  <- exp(tka + eta.ka)
      cl  <- exp(tcl + eta.cl)
      v   <- exp(tv  + eta.v)
      lag(depot) <- tlag
      linCmt() ~ add(add.sd)
    })
  }

  saemControlFast <- saemControl(print = 0, nBurn = 5, nEm = 5, nmc = 1,
                                  nu = c(2, 2, 2))

  test_that("SAEM with logit-bounded param keeps estimate within bounds (#496)", {
    fit <- suppressMessages(suppressWarnings(
      nlmixr(.logitModel, theo_sd, est = "saem", control = saemControlFast)
    ))
    # td1 should be in [0, 1] and named "td1" (not "td1_untransformed")
    expect_true("td1" %in% names(fit$theta))
    expect_false(any(grepl("_untransformed$", names(fit$theta))))
    td1Val <- fit$theta["td1"]
    expect_true(td1Val >= 0 && td1Val <= 1)
  })

  test_that("SAEM with lower-bound-only param (tlag >= 0) works", {
    fit <- suppressMessages(suppressWarnings(
      nlmixr(.lowerModel, theo_sd, est = "saem", control = saemControlFast)
    ))
    expect_true("tlag" %in% names(fit$theta))
    expect_false(any(grepl("_untransformed$", names(fit$theta))))
    tlagVal <- fit$theta["tlag"]
    expect_true(tlagVal >= 0)
  })

  test_that("unbounded params are not transformed (regression)", {
    .unboundedModel <- function() {
      ini({
        tka   <- 0.45
        tcl   <- 1.0
        tv    <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v  ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv  + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    fit <- suppressMessages(suppressWarnings(
      nlmixr(.unboundedModel, theo_sd, est = "saem", control = saemControlFast)
    ))
    expect_false(any(grepl("_untransformed$", names(fit$theta))))
    expect_true(all(c("tka", "tcl", "tv") %in% names(fit$theta)))
  })

  test_that("fit$ui is restored to original model after SAEM", {
    fit <- suppressMessages(suppressWarnings(
      nlmixr(.logitModel, theo_sd, est = "saem", control = saemControlFast)
    ))
    .iniDf <- fit$ui$iniDf
    .thetaRows <- .iniDf[is.na(.iniDf$neta1), ]
    # Original param name restored
    expect_true("td1" %in% .thetaRows$name)
    expect_false(any(grepl("_untransformed$", .thetaRows$name)))
    # Original bounds restored
    .td1 <- .thetaRows[.thetaRows$name == "td1", ]
    expect_equal(.td1$lower, 0)
    expect_equal(.td1$upper, 1)
  })

  test_that("control$boundedTransform = FALSE disables the transform", {
    fit <- suppressMessages(suppressWarnings(
      nlmixr(.logitModel, theo_sd, est = "saem",
             control = saemControl(print = 0, nBurn = 5, nEm = 5,
                                    nmc = 1, nu = c(2, 2, 2),
                                    boundedTransform = FALSE))
    ))
    # Without the transform, td1 should still be named td1
    # (no rewriting happened) but may go out of bounds
    expect_true("td1" %in% names(fit$theta))
    expect_false(any(grepl("_untransformed$", names(fit$theta))))
  })

  test_that("fixed params are not transformed", {
    .m <- function() {
      ini({
        tka   <- 0.45
        tcl   <- 1.0
        tv    <- 3.45
        td1   <- fix(0.5)
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v  ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv  + eta.v)
        d1 <- td1
        f(depot) <- d1
        linCmt() ~ add(add.sd)
      })
    }
    fit <- suppressMessages(suppressWarnings(
      nlmixr(.m, theo_sd, est = "saem", control = saemControlFast)
    ))
    expect_true("td1" %in% names(fit$theta))
    expect_equal(setNames(fit$theta["td1"], NULL), 0.5)
  })

  foceiControlFast <- foceiControl(print = 0, maxOuterIterations = 0L)

  test_that("FOCEI with logit-bounded param keeps estimate within bounds", {
    fit <- suppressMessages(suppressWarnings(
      nlmixr(.logitModel, theo_sd, est = "focei", control = foceiControlFast)
    ))
    expect_true("td1" %in% names(fit$theta))
    expect_false(any(grepl("_untransformed$", names(fit$theta))))
    td1Val <- fit$theta["td1"]
    expect_true(td1Val >= 0 && td1Val <= 1)
  })

  test_that("FOCEI with boundedTransform = FALSE disables the transform", {
    fit <- suppressMessages(suppressWarnings(
      nlmixr(.logitModel, theo_sd, est = "focei",
             control = foceiControl(print = 0, maxOuterIterations = 0L,
                                     boundedTransform = FALSE))
    ))
    expect_true("td1" %in% names(fit$theta))
    expect_false(any(grepl("_untransformed$", names(fit$theta))))
  })

  test_that("boundedTransform arg wired into all 8 control functions", {
    expect_true(saemControl(boundedTransform = FALSE)$boundedTransform == FALSE)
    expect_true(foceiControl(boundedTransform = FALSE)$boundedTransform == FALSE)
    expect_true(n1qn1Control(boundedTransform = FALSE)$boundedTransform == FALSE)
    expect_true(uobyqaControl(boundedTransform = FALSE)$boundedTransform == FALSE)
    expect_true(newuoaControl(boundedTransform = FALSE)$boundedTransform == FALSE)
    expect_true(nlmControl(boundedTransform = FALSE)$boundedTransform == FALSE)
    expect_true(nlsControl(boundedTransform = FALSE)$boundedTransform == FALSE)
    expect_true(optimControl(boundedTransform = FALSE)$boundedTransform == FALSE)
    expect_true(saemControl()$boundedTransform)
    expect_true(foceiControl()$boundedTransform)
    expect_true(n1qn1Control()$boundedTransform)
  })

  rxode2::rxUnloadAll()
})
