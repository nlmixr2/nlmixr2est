nmTest({
  # Weekly-batch sweep: the decoupled "sa"/"imp" covariances applied at fit time
  # (covMethod=) across families and post-hoc via setCov().  Multiple saem/imp
  # fits -> slow; kept out of the essential push/PR subset via .slowBatches.
  skip_on_cran()

  .lc <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  .d <- nlmixr2data::theo_sd

  .isPdFinite <- function(m) {
    is.matrix(m) &&
      all(is.finite(m)) &&
      all(diag(m) > 0) &&
      min(eigen(m, symmetric = TRUE, only.values = TRUE)$values) > 0
  }

  test_that("focei covMethod='sa'/'imp' installs the decoupled covariance", {
    .fs <- suppressWarnings(nlmixr2(.lc, .d, est = "focei", control = foceiControl(print = 0L, covMethod = "sa")))
    expect_equal(.fs$covMethod, "sa")
    expect_true(.isPdFinite(.fs$cov))
    expect_true(all(is.finite(.fs$parFixedDf[["SE"]])))

    .fi <- suppressWarnings(nlmixr2(.lc, .d, est = "focei", control = foceiControl(print = 0L, covMethod = "imp")))
    expect_equal(.fi$covMethod, "imp")
    expect_true(.isPdFinite(.fi$cov))
  })

  # the formatted $parFixed must track the installed cov (issue #816: only
  # $parFixedDf was refreshed, leaving the displayed SE stale)
  .expectParFixedTracksCov <- function(.f, .n = "add.sd") {
    expect_equal(unname(.f$parFixedDf[.n, "SE"]), unname(sqrt(diag(.f$cov))[.n]))
    .seNum <- suppressWarnings(as.numeric(.f$parFixed[.n, "SE"]))
    expect_true(is.finite(.seNum))
    expect_equal(.seNum, signif(unname(.f$parFixedDf[.n, "SE"]), 3), tolerance = 1e-2)
  }

  test_that("setCov() switches any completed fit to sa/imp", {
    .f <- suppressWarnings(nlmixr2(.lc, .d, est = "focei", control = foceiControl(print = 0L)))
    expect_equal(.f$covMethod, "r,s (full)")
    .cov0 <- .f$cov
    suppressMessages(setCov(.f, "sa"))
    expect_equal(.f$covMethod, "sa")
    expect_true(.isPdFinite(.f$cov))
    .expectParFixedTracksCov(.f)
    suppressMessages(setCov(.f, "imp"))
    expect_equal(.f$covMethod, "imp")
    .expectParFixedTracksCov(.f)
    # the original r,s covariance stays recoverable from the cache
    suppressMessages(setCov(.f, "r,s (full)"))
    expect_equal(.f$covMethod, "r,s (full)")
    expect_equal(unname(.f$cov), unname(.cov0))
    # ... and so does the theta-only shape the fit cached alongside it
    suppressMessages(setCov(.f, "r,s"))
    expect_equal(.f$covMethod, "r,s")
  })

  test_that("saem accepts the foreign imp covariance", {
    .f <- suppressWarnings(nlmixr2(
      .lc,
      .d,
      est = "saem",
      control = saemControl(print = 0L, nBurn = 100, nEm = 100, covMethod = "imp")
    ))
    expect_equal(.f$covMethod, "imp")
    expect_true(.isPdFinite(.f$cov))
  })

  test_that("the sa and imp recomputes run at the fit's estimates without moving them", {
    # The "sa" engine's nBurn + nEm SAEM iterations and the "imp" engine's EM step
    # start from the pinned estimates and must leave every parameter there, so
    # the covariance is the one at those estimates.
    .f <- suppressWarnings(nlmixr2(.lc, .d, est = "focei", control = foceiControl(print = 0L, covMethod = "")))
    .a <- .covPinnedRefitArgs(.f)
    .eta <- c("eta.ka", "eta.cl", "eta.v")
    .sa <- suppressWarnings(suppressMessages(nlmixr2(
      .a$ui,
      .a$data,
      est = "saem",
      control = .covEngineControl("sa", saControl(nBurn = 20L, nEm = 20L, nSaCov = 50L))
    )))
    expect_identical(.sa$covMethod, "sa")
    # the hold is internal to the run, never part of the stored control
    expect_false("saemHoldPar" %in% names(.sa$control))
    expect_false("saemHoldPar" %in% names(.sa$env$saemControl))
    expect_equal(.sa$theta[names(.f$theta)], .f$theta, tolerance = 1e-10)
    expect_equal(.sa$omega[.eta, .eta], .f$omega[.eta, .eta], tolerance = 1e-10)
    .ctl <- .covEngineControl("imp", impCovControl(nIter = 3L, isample = 100L))
    .ctl$etaMat <- .a$etaMat
    .imp <- suppressWarnings(suppressMessages(nlmixr2(.a$ui, .a$data, est = "imp", control = .ctl)))
    expect_true(is.matrix(.imp$env$impCovInternal))
    expect_false("impFrozen" %in% names(.imp$control))
    expect_false("impFrozen" %in% names(.imp$env$impmapControl))
    expect_equal(.imp$theta[names(.f$theta)], .f$theta, tolerance = 1e-10)
    expect_equal(.imp$omega[.eta, .eta], .f$omega[.eta, .eta], tolerance = 1e-10)
  })

  test_that("setCov(\"sa\") on a SAEM fit continues its chains with no warm-up iterations", {
    .f <- suppressWarnings(nlmixr2(.lc, .d, est = "saem", control = saemControl(print = 0L, nBurn = 100, nEm = 100)))
    expect_identical(.f$covMethod, "sa")
    .native <- sqrt(diag(.f$cov))
    .all <- .f$phiM
    expect_identical(dim(.all), c(12L, 3L, 200L, 3L))
    .st <- .saemChainState(.f)
    expect_identical(dim(.st$phiM), c(36L, 3L))
    # row i + k * N is subject i of chain k
    expect_identical(.st$phiM[2L + 12L, ], .all[2L, 2L, 200L, ])
    # the sigma2 the fit's Louis residual score read: add.sd^2 at convergence
    expect_equal(.st$sigma2, unname(.f$theta["add.sd"])^2, tolerance = 1e-3)
    expect_identical(.st$mpostPhi, .f$saem$mpost_phi)
    # a fit without chains starts cold
    .fo <- suppressWarnings(nlmixr2(.lc, .d, est = "focei", control = foceiControl(print = 0L, covMethod = "")))
    expect_null(.saemChainState(.fo))

    .eta <- c("eta.ka", "eta.cl", "eta.v")
    .a <- .covPinnedRefitArgs(.f)
    .sa <- suppressWarnings(suppressMessages(nlmixr2(
      .a$ui,
      .a$data,
      est = "saem",
      control = .covEngineControl("sa", saControl(), .st)
    )))
    expect_identical(.sa$covMethod, "sa")
    expect_false("saemWarmState" %in% names(.sa$control))
    # no estimation iteration ran, so nothing moved
    expect_identical(nrow(.sa$parHistData), 0L)
    expect_equal(.sa$theta[names(.f$theta)], .f$theta, tolerance = 1e-10)
    expect_equal(.sa$omega[.eta, .eta], .f$omega[.eta, .eta], tolerance = 1e-10)
    # the same estimand as the fit's own sa phase, to Monte Carlo noise (about
    # 1% here); the residual SE depends on the seeded residual statistic, which
    # left at the kernel's start gave an SE of about 9
    expect_equal(sqrt(diag(.sa$cov)), .native, tolerance = 0.1)

    suppressMessages(setCov(.f, "sa", saControl(nSaCov = 400L)))
    expect_identical(.f$covMethod, "sa")
    expect_equal(sqrt(diag(.f$cov)), .native, tolerance = 0.1)
  })

  test_that("impmap accepts the foreign sa covariance", {
    .f <- suppressWarnings(nlmixr2(.lc, .d, est = "impmap", control = impmapControl(print = 0L, covMethod = "sa")))
    expect_equal(.f$covMethod, "sa")
    expect_true(.isPdFinite(.f$cov))
  })
})
