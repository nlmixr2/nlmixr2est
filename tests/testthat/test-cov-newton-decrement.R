test_that("the decrement is 0.5 g' H^-1 g of a quadratic OFV", {
  # OFV(x) = 0.5 (x - m)' H (x - m); the covariance of -2 log L is 2 H^-1
  .h <- matrix(c(4, 1, 0, 1, 3, 0.5, 0, 0.5, 2), 3, 3)
  .nm <- c("a", "b", "c")
  dimnames(.h) <- list(.nm, .nm)
  .m <- c(a = 1, b = -2, c = 0.5)
  .x <- c(a = 1.4, b = -1.7, c = 0.2)
  .g <- drop(.h %*% (.x - .m))
  names(.g) <- .nm
  .cov <- 2 * solve(.h)
  .ofv <- 0.5 * drop(crossprod(.x - .m, .h %*% (.x - .m)))
  expect_equal(.covNewtonDecrement(.g, .cov), .ofv)
  expect_equal(.covNewtonDecrement(.g, .cov), 0.5 * drop(crossprod(.g, solve(.h, .g))))
  # parameters missing from either side, or not finite, are left out
  .g2 <- c(.g, d = 3)
  expect_equal(.covNewtonDecrement(.g2, .cov), .ofv)
  .cov2 <- .cov
  .cov2["c", ] <- .cov2[, "c"] <- NA_real_
  expect_equal(.covNewtonDecrement(.g, .cov2), 0.25 * drop(crossprod(.g[1:2], .cov[1:2, 1:2] %*% .g[1:2])))
  expect_identical(.covNewtonDecrement(NULL, .cov), NA_real_)
  expect_identical(.covNewtonDecrement(c(z = 1), .cov), NA_real_)
})

test_that("an R-only covariance stands in for a sandwich", {
  .e <- new.env(parent = emptyenv())
  .rs <- matrix(2, 1, 1, dimnames = list("a", "a"))
  .fit <- list(cov = .rs)
  expect_identical(.covNewtonCov(.fit, .e), .rs)
  .e$covR <- matrix(1, 1, 1, dimnames = list("a", "a"))
  expect_identical(.covNewtonCov(.fit, .e), .e$covR)
  .s <- new.env(parent = emptyenv())
  expect_null(.covNewtonCov(list(cov = .rs, covMethod = "s"), .s))
  expect_null(.covNewtonCov(list(cov = .rs, covMethod = "|s|"), .s))
  expect_identical(.covNewtonCov(list(cov = .rs, covMethod = "sa"), .s), .rs)
  expect_identical(.covNewtonCov(list(cov = .rs, covMethod = "linFim"), .s), .rs)
  .s$covList <- list(r = .e$covR)
  expect_identical(.covNewtonCov(list(cov = .rs, covMethod = "r,s"), .s), .e$covR)
  .s$covList <- list(`r (full)` = .e$covR)
  expect_identical(.covNewtonCov(list(cov = .rs, covMethod = "r,s (full)"), .s), .e$covR)
})

test_that("FOCEi-family methods are not checked", {
  expect_false(.covNewtonEst("focei"))
  expect_false(.covNewtonEst("laplace"))
  expect_false(.covNewtonEst("posthoc"))
  expect_true(.covNewtonEst("saem"))
  expect_true(.covNewtonEst("nlme"))
  expect_false(.covNewtonEst(""))
})

test_that("the warning fits on one runInfo line", {
  expect_lt(nchar(.covNewtonWarnMsg(123456789)), 75)
  expect_identical(
    .covNewtonWarnMsg(4.2134),
    "SEs at a non-stationary estimate (decrement 4.21); polish with focei"
  )
})

test_that("covNewtonWarn is validated", {
  expect_equal(saemControl()$covNewtonWarn, 1)
  expect_true(is.na(saemControl(covNewtonWarn = NA)$covNewtonWarn))
  expect_error(saemControl(covNewtonWarn = -1))
})

nmTest({
  test_that("a truncated SAEM fit warns; a converged one and FOCEi do not", {
    skip_on_cran()
    one.cmt <- function() {
      ini({
        tka <- 1.2
        tcl <- 2
        tv <- 2.5
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
    .d <- nlmixr2data::theo_sd
    .saem <- function(nBurn, nEm, warn = 1) {
      suppressMessages(nlmixr2(
        one.cmt,
        .d,
        "saem",
        saemControl(print = 0, nBurn = nBurn, nEm = nEm, covNewtonWarn = warn)
      ))
    }
    .re <- "^SEs at a non-stationary estimate"
    .short <- .saem(3, 2)
    expect_gt(.short$covNewtonDecrement, 1)
    expect_length(grep(.re, .short$runInfo), 1L)
    .off <- .saem(3, 2, NA)
    expect_null(.off$covNewtonDecrement)
    expect_false(any(grepl(.re, .off$runInfo)))
    .conv <- .saem(200, 100)
    expect_lt(.conv$covNewtonDecrement, 1)
    expect_false(any(grepl(.re, .conv$runInfo)))
    .focei <- suppressMessages(nlmixr2(.conv, .d, "focei", foceiControl(print = 0)))
    expect_null(.focei$covNewtonDecrement)
    expect_false(any(grepl(.re, .focei$runInfo)))
  })
})

nmTest({
  test_that("models outside the analytic gradient's scope use central differences", {
    skip_on_cran()
    .ode <- function() {
      ini({
        tka <- 1.2
        tcl <- 2
        tv <- 2.5
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
    .lin <- function() {
      ini({
        tka <- 1.2
        tcl <- 2
        tv <- 2.5
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
    .d <- nlmixr2data::theo_sd
    .saem <- function(mod, warn = 1) {
      suppressMessages(nlmixr2(mod, .d, "saem", saemControl(print = 0, nBurn = 3, nEm = 2, covNewtonWarn = warn)))
    }
    .f <- .saem(.ode, NA)
    .an <- .covNewtonGrad(.f, .f$cov)
    .fd <- .covNewtonGradFd(.f, .f$cov)
    # only the analytic route has the Omega (Cholesky) elements
    expect_true(any(grepl("^om[.]chol", names(.an))))
    expect_named(.fd, c("tka", "tcl", "tv", "add.sd"))
    expect_equal(unname(.fd), unname(.an[names(.fd)]), tolerance = 0.02)
    .l <- .saem(.lin)
    expect_gt(.l$covNewtonDecrement, 1)
    expect_length(grep("^SEs at a non-stationary estimate", .l$runInfo), 1L)
  })
})

test_that("a rejected probe falls back to a one-sided difference", {
  .ui <- rxode2::rxUiDecompress(rxode2::rxode2(function() {
    ini({
      p1 <- 0.9995
      tcl <- 1
      eta.cl ~ 0.1
      add.sd <- 0.7
    })
    model({
      cl <- exp(tcl + eta.cl) * p1
      cl ~ add(add.sd)
    })
  }))
  .fit <- list(finalUi = .ui)
  .cov <- diag(2)
  dimnames(.cov) <- list(c("p1", "tcl"), c("p1", "tcl"))
  # a quadratic OFV in place of the nested FOCEi pass; p1 above 1 is rejected
  local_mocked_bindings(
    .covNewtonCtl = function(fit, fast) NULL,
    .covNewtonFocei = function(ui, fit, ctl) {
      .e <- stats::setNames(ui$iniDf$est, ui$iniDf$name)
      if (.e[["p1"]] > 1) {
        stop("probability above 1")
      }
      list(objf = 3 * .e[["p1"]]^2 + 2 * .e[["tcl"]]^2)
    }
  )
  .g <- .covNewtonGradFd(.fit, .cov)
  expect_equal(.g[["tcl"]], 4, tolerance = 1e-6)
  expect_equal(.g[["p1"]], 6 * 0.9995, tolerance = 1e-2)
  # both probes rejected at the default step: a smaller one
  local_mocked_bindings(
    .covNewtonFocei = function(ui, fit, ctl) {
      .e <- stats::setNames(ui$iniDf$est, ui$iniDf$name)
      if (abs(.e[["p1"]] - 0.9995) > 5e-4) {
        stop("out of range")
      }
      list(objf = 3 * .e[["p1"]]^2 + 2 * .e[["tcl"]]^2)
    }
  )
  .g <- .covNewtonGradFd(.fit, .cov)
  expect_equal(.g[["p1"]], 6 * 0.9995, tolerance = 1e-6)
  # and no partial gradient when a direction cannot be differenced at all
  local_mocked_bindings(
    .covNewtonFocei = function(ui, fit, ctl) {
      .e <- stats::setNames(ui$iniDf$est, ui$iniDf$name)
      if (.e[["p1"]] != 0.9995) {
        stop("out of range")
      }
      list(objf = 3 * .e[["p1"]]^2 + 2 * .e[["tcl"]]^2)
    }
  )
  expect_null(.covNewtonGradFd(.fit, .cov))
})
