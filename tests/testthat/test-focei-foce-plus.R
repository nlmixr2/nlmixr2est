nmTest({
  # foce = "foce+" revives the pre-6.1.0 FOCE that keeps the live conditional
  # residual variance R (vs. the default "nonmem" FOCE that freezes R at eta=0).
  # For a proportional/combined-error model the two objectives differ; for pure
  # additive error R does not depend on eta, so they must coincide.
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")

  one.cmt <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7; prop.sd <- 0.1
          eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd) + prop(prop.sd) })
  }

  add.cmt <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
          eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd) })
  }

  d <- nlmixr2data::theo_sd

  test_that("foce=\"foce+\" runs and differs from \"nonmem\" for proportional error", {
    fn <- suppressWarnings(suppressMessages(
      nlmixr(one.cmt, d, "foce",
             foceiControl(print = 0L, calcTables = FALSE, covMethod = ""))))
    fp <- suppressWarnings(suppressMessages(
      nlmixr(one.cmt, d, "foce",
             foceiControl(print = 0L, calcTables = FALSE, covMethod = "", foce = "foce+"))))
    expect_true(is.finite(fn$objective))
    expect_true(is.finite(fp$objective))
    expect_false(isTRUE(all.equal(fn$objective, fp$objective, tolerance = 1e-4)))
  })

  test_that("foce=\"foce+\" equals \"nonmem\" for pure additive error", {
    fn <- suppressWarnings(suppressMessages(
      nlmixr(add.cmt, d, "foce",
             foceiControl(print = 0L, calcTables = FALSE, covMethod = ""))))
    fp <- suppressWarnings(suppressMessages(
      nlmixr(add.cmt, d, "foce",
             foceiControl(print = 0L, calcTables = FALSE, covMethod = "", foce = "foce+"))))
    expect_equal(fn$objective, fp$objective, tolerance = 1e-4)
  })

  test_that("focep/mfocep/ifocep equal foce with foce=\"foce+\"", {
    ref <- suppressWarnings(suppressMessages(
      nlmixr(one.cmt, d, "foce",
             foceiControl(print = 0L, calcTables = FALSE, covMethod = "", foce = "foce+"))))
    for (est in c("focep", "mfocep", "ifocep")) {
      fit <- suppressWarnings(suppressMessages(
        nlmixr(one.cmt, d, est,
               foceiControl(print = 0L, calcTables = FALSE, covMethod = ""))))
      expect_true(is.finite(fit$objective))
      expect_equal(fit$objective, ref$objective, tolerance = 1e-3)
    }
  })
})
