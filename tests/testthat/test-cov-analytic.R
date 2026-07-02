# Analytic FOCEI/FOCE covariance (foceiCov / foceiCovAnalytic) via the uniform
# direction assembly: every mu-ref theta reuses its eta's direction, every other
# structural theta (covariate coefficient OR eta-less) gets its own true-sensitivity
# direction.  covFull=FALSE (default) is the theta-only cov; covFull=TRUE spans
# theta + sigma + Omega.  Runs on a released rxode2 (needs rxExpandSens2_/symengine).

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

test_that("foceiCov: default is theta-only, covFull=TRUE spans theta+sigma+Omega", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  fit <- suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "")))
  # default covFull=FALSE -> structural thetas only (the original nlmixr2est scope)
  r0 <- foceiCov(fit)
  expect_false(is.null(r0))
  expect_identical(r0$method, "analytic")
  expect_setequal(r0$params, c("tka", "tcl", "tv"))
  # covFull=TRUE -> full theta + sigma + Omega
  r <- foceiCov(fit, covFull = TRUE)
  expect_setequal(r$params, c("tka", "tcl", "tv", "add.sd", "om.tka", "om.tcl", "om.tv"))
  expect_true(all(is.finite(r$se)) && all(r$se > 0))
  # structural-theta SEs match NONMEM $COV MATRIX=R (0.19180/0.08352/0.04661)
  expect_equal(unname(r$se[c("tka", "tcl", "tv")]),
               c(0.18868, 0.08351, 0.04617), tolerance = 0.03)
  expect_lt(max(abs(r$R - t(r$R))), 1e-6)                        # R symmetric
  expect_true(all(eigen(r$R, symmetric = TRUE, only.values = TRUE)$values > 0))  # PD
})

test_that("the fd2 tier reproduces the exact 3rd-order tier (full covariance)", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  fit <- suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "")))
  r3 <- foceiCovAnalytic(fit, sens = "exact3", covFull = TRUE)
  r2 <- foceiCovAnalytic(fit, sens = "fd2",    covFull = TRUE)
  expect_false(is.null(r2))
  expect_setequal(r2$params, r3$params)
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
  r <- foceiCov(fit, covFull = TRUE)
  expect_false(is.null(r))
  expect_identical(r$method, "analytic")             # block Omega via the E-basis derivatives
  expect_true(any(grepl("^cov\\.", r$params)))       # an off-diagonal Omega covariance SE
  expect_true(all(is.finite(r$se)) && all(r$se > 0))
  expect_true(all(eigen(r$cov, symmetric = TRUE, only.values = TRUE)$values > 0))
  r2 <- foceiCovAnalytic(fit, sens = "fd2", covFull = TRUE)
  expect_equal(unname(r$se), unname(r2$se[r$params]), tolerance = 1e-3)
})

test_that("single random-effect model is handled analytically (no sapply collapse)", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  one.eta <- function() {                            # exactly one eta, all thetas mu-referenced
    ini({ tcl <- log(2.7); eta.cl ~ 0.1; add.sd <- 0.7 })
    model({
      ka <- 1.5; cl <- exp(tcl + eta.cl); v <- 31.5
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  fit <- suppressMessages(nlmixr(one.eta, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "")))
  r <- foceiCovAnalytic(fit, covFull = TRUE)
  expect_false(is.null(r))
  expect_setequal(r$params, c("tcl", "add.sd", "om.tcl"))
  expect_true(all(is.finite(r$se)) && all(r$se > 0))
})

test_that("covariate coefficient handled via its own direction (matches FD)", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  covm <- function() {                               # wt_eff is a covariate coefficient on cl
    ini({
      tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5); wt_eff <- 0.75
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + wt_eff * log(WT / 70) + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  # analytic on a plain fit
  fa <- suppressMessages(nlmixr(covm, nlmixr2data::theo_sd, "focei",
                                foceiControl(print = 0L, covMethod = "", sigdig = 6)))
  r <- foceiCovAnalytic(fa, covFull = TRUE)
  expect_false(is.null(r))
  expect_true("wt_eff" %in% r$params)                # covariate coefficient present
  expect_true(all(is.finite(r$se)) && all(r$se > 0))
  # value-check against the finite-difference full covariance
  ffd <- suppressMessages(nlmixr(covm, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "r", covFull = TRUE, sigdig = 6)))
  seFD <- sqrt(abs(diag(ffd$cov)))
  expect_equal(unname(r$se[names(seFD)]), unname(seFD), tolerance = 5e-3)
})

test_that("eta-less structural theta handled via its own direction (matches FD)", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  elm <- function() {                                # tv has NO eta.v -> eta-less structural theta
    ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5); eta.ka ~ 0.6; eta.cl ~ 0.3; add.sd <- 0.7 })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  fa <- suppressMessages(nlmixr(elm, nlmixr2data::theo_sd, "focei",
                                foceiControl(print = 0L, covMethod = "", sigdig = 6)))
  r <- foceiCovAnalytic(fa, covFull = TRUE)
  expect_false(is.null(r))
  expect_true("tv" %in% r$params)                    # eta-less theta present (would bow out under scaling)
  expect_true(all(is.finite(r$se)) && all(r$se > 0))
  ffd <- suppressMessages(nlmixr(elm, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "r", covFull = TRUE, sigdig = 6)))
  seFD <- sqrt(abs(diag(ffd$cov)))
  expect_equal(unname(r$se[names(seFD)]), unname(seFD), tolerance = 5e-3)
})

test_that("covType='analytic' installs the analytic cov end-to-end", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  fa <- suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
                                foceiControl(print = 0L, covMethod = "r", covType = "analytic",
                                             covFull = TRUE, sigdig = 6)))
  # $cov spans theta+sigma+Omega and matches the FD engine within tolerance
  expect_setequal(colnames(fa$cov), c("tka", "tcl", "tv", "add.sd", "om.tka", "om.tcl", "om.tv"))
  ffd <- suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "r", covType = "fd",
                                              covFull = TRUE, sigdig = 6)))
  seA <- sqrt(abs(diag(fa$cov))); seF <- sqrt(abs(diag(ffd$cov)))
  expect_equal(unname(seA[names(seF)]), unname(seF), tolerance = 5e-3)
})
