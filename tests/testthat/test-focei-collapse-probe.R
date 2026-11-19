test_that("collapseProbe is validated and off by default", {
  expect_false(foceiControl()$collapseProbe)
  expect_true(foceiControl(collapseProbe = TRUE)$collapseProbe)
  expect_error(foceiControl(collapseProbe = NA))
  expect_error(foceiControl(collapseProbe = "yes"))
})

test_that("the probe omega raises collapsed etas and drops their covariances", {
  .om <- matrix(c(1, 0.01, 0.01, 1e-6), 2, 2, dimnames = list(c("eta.a", "eta.b"), c("eta.a", "eta.b")))
  .p <- .foceiCollapseProbeOmega(.om, c(eta.a = 2, eta.b = 300), "eta.b")
  expect_equal(.p, matrix(c(1, 0, 0, 300), 2, 2, dimnames = dimnames(.om)))
})

test_that("no probe snapshot unless requested on an estimating mixed fit", {
  .env <- new.env(parent = emptyenv())
  .env$rxInv <- structure(new.env(parent = emptyenv()), class = "rxSymInvCholEnv")
  .env$thetaIni <- c(`THETA[1]` = 1)
  expect_null(.foceiCollapseProbeSnapshot(.env, list(collapseProbe = FALSE, maxOuterIterations = 10L)))
  expect_null(.foceiCollapseProbeSnapshot(.env, list(collapseProbe = TRUE, maxOuterIterations = 0L)))
  expect_named(
    .foceiCollapseProbeSnapshot(.env, list(collapseProbe = TRUE, maxOuterIterations = 10L)),
    c("rxInv", "thetaIni"),
    ignore.order = TRUE
  )
  .env$mixIdx <- 1L
  expect_null(.foceiCollapseProbeSnapshot(.env, list(collapseProbe = TRUE, maxOuterIterations = 10L)))
  # without a snapshot, or without a collapsed omega, the fit is returned as is
  .fit <- new.env(parent = emptyenv())
  .fit$omega <- matrix(1, 1, 1, dimnames = list("eta.a", "eta.a"))
  expect_identical(.foceiCollapseProbe(.fit, NULL, NULL, list()), .fit)
  .ui <- list(iniDf = data.frame(name = "eta.a", neta1 = 1, neta2 = 1, est = 1, fix = FALSE))
  expect_identical(.foceiCollapseProbe(.fit, list(thetaIni = 1), .ui, list()), .fit)
})

test_that("collapseProbe re-fits a collapsed random-slope omega (#1394)", {
  skip_on_cran()
  set.seed(1394)
  .d <- do.call(
    rbind,
    lapply(seq_len(30), function(i) {
      .b0 <- 250 + rnorm(1, 0, sqrt(600))
      .b1 <- 10 + rnorm(1, 0, sqrt(35))
      data.frame(ID = i, TIME = 0:9, DV = .b0 + .b1 * (0:9) + rnorm(10, 0, 25))
    })
  )
  .mod <- function() {
    ini({
      tint <- 200
      tsl <- 5
      eta.int ~ 400
      eta.sl ~ 300
      add.sd <- 100
    })
    model({
      b0 <- tint + eta.int
      b1 <- tsl + eta.sl
      y <- b0 + b1 * time
      y ~ add(add.sd)
    })
  }
  .ctl <- function(...) foceiControl(print = 0, covMethod = "", calcTables = FALSE, ...)
  # the interior optimum, from the default omega bound
  .ref <- suppressMessages(nlmixr2(.mod, .d, "focei", .ctl()))
  # without the bound the slope variance collapses toward 0
  .collapsed <- suppressMessages(nlmixr2(.mod, .d, "focei", .ctl(diagOmegaBoundLower = 1e8)))
  expect_lt(.collapsed$omega["eta.sl", "eta.sl"], 0.3)
  expect_gt(.collapsed$objective, .ref$objective + 50)
  .probed <- suppressMessages(nlmixr2(
    .mod,
    .d,
    "focei",
    foceiControl(print = 0, covMethod = "r", diagOmegaBoundLower = 1e8, collapseProbe = TRUE)
  ))
  expect_lt(abs(.probed$objective - .ref$objective), 1)
  expect_gt(.probed$omega["eta.sl", "eta.sl"], 10)
  expect_true(any(grepl("^omega collapse probe lowered objective by [0-9.]+: eta.sl$", .probed$runInfo)))
  expect_false(any(grepl("^omega collapsed", .probed$runInfo)))
  # the kept fit is finished from the re-fit: tables and covariance included
  expect_equal(nrow(.probed), nrow(.d))
  expect_true(all(is.finite(.probed$CWRES)))
  expect_true(is.matrix(.probed$cov))
})
