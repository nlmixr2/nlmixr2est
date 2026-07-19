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
    is.matrix(m) && all(is.finite(m)) && all(diag(m) > 0) &&
      min(eigen(m, symmetric = TRUE, only.values = TRUE)$values) > 0
  }

  test_that("focei covMethod='sa'/'imp' installs the decoupled covariance", {
    .fs <- suppressWarnings(nlmixr2(.lc, .d, est = "focei",
                                    control = foceiControl(print = 0L, covMethod = "sa")))
    expect_equal(.fs$covMethod, "sa")
    expect_true(.isPdFinite(.fs$cov))
    expect_true(all(is.finite(.fs$parFixedDf[["SE"]])))

    .fi <- suppressWarnings(nlmixr2(.lc, .d, est = "focei",
                                    control = foceiControl(print = 0L, covMethod = "imp")))
    expect_equal(.fi$covMethod, "imp")
    expect_true(.isPdFinite(.fi$cov))
  })

  test_that("setCov() switches any completed fit to sa/imp", {
    .f <- suppressWarnings(nlmixr2(.lc, .d, est = "focei",
                                   control = foceiControl(print = 0L)))
    expect_equal(.f$covMethod, "r,s")
    suppressMessages(setCov(.f, "sa"))
    expect_equal(.f$covMethod, "sa")
    expect_true(.isPdFinite(.f$cov))
    suppressMessages(setCov(.f, "imp"))
    expect_equal(.f$covMethod, "imp")
    # the original r,s covariance stays recoverable from the cache
    suppressMessages(setCov(.f, "r,s"))
    expect_equal(.f$covMethod, "r,s")
  })

  test_that("saem accepts the foreign imp covariance", {
    .f <- suppressWarnings(nlmixr2(.lc, .d, est = "saem",
                                   control = saemControl(print = 0L, nBurn = 100, nEm = 100,
                                                         covMethod = "imp")))
    expect_equal(.f$covMethod, "imp")
    expect_true(.isPdFinite(.f$cov))
  })

  test_that("impmap accepts the foreign sa covariance", {
    .f <- suppressWarnings(nlmixr2(.lc, .d, est = "impmap",
                                   control = impmapControl(print = 0L, covMethod = "sa")))
    expect_equal(.f$covMethod, "sa")
    expect_true(.isPdFinite(.f$cov))
  })
})
