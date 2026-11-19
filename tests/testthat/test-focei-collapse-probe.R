test_that("collapseProbe is validated and off by default", {
  expect_false(foceiControl()$collapseProbe)
  expect_true(foceiControl(collapseProbe = TRUE)$collapseProbe)
  expect_error(foceiControl(collapseProbe = NA))
  expect_error(foceiControl(collapseProbe = "yes"))
})

test_that("the probe omega raises collapsed etas and keeps the block", {
  .nm <- c("eta.a", "eta.b", "eta.c")
  .om <- matrix(c(1, 0.008, 0, 0.008, 1e-4, 0, 0, 0, 2), 3, 3, dimnames = list(.nm, .nm))
  .p <- .foceiCollapseProbeOmega(.om, c(eta.a = 2, eta.b = 300, eta.c = 1), "eta.b")
  expect_equal(diag(.p), c(eta.a = 1, eta.b = 300, eta.c = 2))
  # the correlation is halved, not dropped, so the omega parameters are kept
  expect_equal(cov2cor(.p)["eta.a", "eta.b"], cov2cor(.om)["eta.a", "eta.b"] / 2)
  expect_identical(.p != 0, .om != 0)
  expect_true(all(eigen(.p, only.values = TRUE)$values > 0))
  expect_identical(
    length(rxode2::rxSymInvCholCreate(.p, diag.xform = "sqrt")$theta),
    length(rxode2::rxSymInvCholCreate(.om, diag.xform = "sqrt")$theta)
  )
  # a covariance the fit left at exactly 0 keeps the starting block
  .om0 <- .om
  .om[1, 2] <- .om[2, 1] <- 0
  expect_identical(.foceiCollapseProbeOmega(.om, c(eta.b = 300), "eta.b", .om0) != 0, .om0 != 0)
})

test_that("the probe start stays inside narrow and infinite bounds", {
  .c <- .foceiCollapseProbeClamp(c(5, 0.0005, -3, 2), c(0, 1e-4, -Inf, 1), c(1, 1e-3, Inf, 3))
  expect_true(all(.c[1:2] > c(0, 1e-4) & .c[1:2] < c(1, 1e-3)))
  expect_identical(.c[3:4], c(-3, 2))
  # a start already inside the bounds is not moved, however close to one
  expect_identical(.foceiCollapseProbeClamp(1e-5, 1e-6, 10), 1e-5)
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
    foceiControl(print = 0, covMethod = "", diagOmegaBoundLower = 1e8, collapseProbe = TRUE)
  ))
  expect_lt(abs(.probed$objective - .ref$objective), 1)
  expect_gt(.probed$omega["eta.sl", "eta.sl"], 10)
  expect_true(any(grepl("^omega collapse probe lowered objective by [0-9.]+: eta.sl$", .probed$runInfo)))
  expect_false(any(grepl("^omega collapsed", .probed$runInfo)))
  # the kept fit is finished from the re-fit, tables included
  expect_equal(nrow(.probed), nrow(.d))
  expect_true(all(is.finite(.probed$CWRES)))
})
