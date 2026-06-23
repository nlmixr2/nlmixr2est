# Analytic FOCEI/FOCE covariance (foceiCovAnalytic / foceiCovFD / foceiCov).
# These run on a released rxode2 (the pure-R 3rd-order sensitivity generator is
# used when rxExpandSens3_ is absent).

test_that("analytic covariance reproduces FOCEI structural SEs on theophylline", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  one.cmt <- function() {
    ini({
      tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  fit <- suppressMessages(nlmixr(one.cmt, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "")))
  r <- foceiCov(fit)
  expect_false(is.null(r))
  expect_identical(r$method, "analytic")
  expect_setequal(r$params, c("tka", "tcl", "tv"))
  # matches NONMEM $COV MATRIX=R (tka 0.19180, tcl 0.08352, tv 0.04661) to ~few %
  expect_equal(unname(r$se[c("tka", "tcl", "tv")]),
               c(0.18868, 0.08351, 0.04617), tolerance = 0.03)
  # R-matrix is symmetric positive-definite
  expect_lt(max(abs(r$R - t(r$R))), 1e-8)
  expect_true(all(eigen(r$R, symmetric = TRUE, only.values = TRUE)$values > 0))
})

test_that("foceiCovFD returns the full theta+Omega+residual covariance", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  one.cmt <- function() {
    ini({
      tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  fit <- suppressMessages(nlmixr(one.cmt, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "")))
  r <- foceiCov(fit, omega = TRUE, residual = TRUE)
  expect_false(is.null(r))
  expect_identical(r$method, "fd")
  # full parameter set: structural thetas + residual + Omega variances
  expect_true(all(c("tka", "tcl", "tv", "add.sd") %in% r$params))
  expect_equal(length(r$params), 7L)              # 3 theta + 1 sigma + 3 omega
  expect_true(all(is.finite(r$se)) && all(r$se > 0))
  # structural-theta SEs agree with the analytic path
  ra <- foceiCov(fit)
  expect_equal(unname(r$se[c("tka", "tcl", "tv")]),
               unname(ra$se[c("tka", "tcl", "tv")]), tolerance = 0.05)
})
