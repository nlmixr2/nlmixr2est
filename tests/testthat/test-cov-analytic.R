# Full analytic FOCEI/FOCE covariance (foceiCov / foceiCovAnalytic) and the
# fallback ladder.  Runs on a released rxode2 (the pure-R 3rd-order sensitivity
# generator is used when rxExpandSens3_ is absent).

.cov_one_cmt <- function() {
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

test_that("foceiCov returns the full theta+sigma+Omega analytic covariance", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  fit <- suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "")))
  r <- foceiCov(fit)
  expect_false(is.null(r))
  expect_identical(r$method, "analytic")
  # every population parameter present: 3 theta + 1 sigma + 3 Omega variances
  expect_setequal(r$params, c("tka", "tcl", "tv", "add.sd", "om.tka", "om.tcl", "om.tv"))
  expect_true(all(is.finite(r$se)) && all(r$se > 0))
  # structural-theta SEs match NONMEM $COV MATRIX=R (0.19180/0.08352/0.04661)
  expect_equal(unname(r$se[c("tka", "tcl", "tv")]),
               c(0.18868, 0.08351, 0.04617), tolerance = 0.03)
  # R matrix symmetric positive-definite
  expect_lt(max(abs(r$R - t(r$R))), 1e-6)
  expect_true(all(eigen(r$R, symmetric = TRUE, only.values = TRUE)$values > 0))
})

test_that("the fd2 tier reproduces the exact 3rd-order tier (full covariance)", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  fit <- suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "")))
  r3 <- foceiCovAnalytic(fit, sens = "exact3")
  r2 <- foceiCovAnalytic(fit, sens = "fd2")
  expect_false(is.null(r2))
  expect_setequal(r2$params, r3$params)
  # Shi finite differences of the analytic 2nd-order sensitivities reproduce the
  # exact 3rd-order R matrix to FD accuracy
  expect_equal(unname(r2$se), unname(r3$se), tolerance = 1e-3)
})

test_that("block Omega is handled analytically, with off-diagonal covariance SEs", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  blk <- function() {
    ini({
      tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
      eta.ka ~ 0.6
      eta.cl + eta.v ~ c(0.3, 0.03, 0.1)            # block between eta.cl and eta.v
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
  fit <- suppressMessages(nlmixr(blk, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "")))
  r <- foceiCov(fit)
  expect_false(is.null(r))
  expect_identical(r$method, "analytic")             # block Omega via the E-basis derivatives
  expect_true(any(grepl("^cov\\.", r$params)))       # an off-diagonal Omega covariance SE
  expect_true(all(is.finite(r$se)) && all(r$se > 0))
  expect_true(all(eigen(r$cov, symmetric = TRUE, only.values = TRUE)$values > 0))
})
