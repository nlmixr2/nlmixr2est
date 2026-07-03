# Analytic FOCEI covariance (.foceiCov / .foceiCovAnalytic) via the uniform
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

test_that(".foceiCov: default is theta-only, covFull=TRUE spans theta+sigma+Omega", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  fit <- suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "")))
  # default covFull=FALSE -> structural thetas only (the original nlmixr2est scope)
  r0 <- .foceiCov(fit)
  expect_false(is.null(r0))
  expect_equal(r0$method, "analytic-fd2")            # fd2 is the default tier
  expect_setequal(r0$params, c("tka", "tcl", "tv"))
  # covFull=TRUE -> full theta + sigma + Omega
  r <- .foceiCov(fit, covFull = TRUE)
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
  r3 <- .foceiCovAnalytic(fit, sens = "exact3", covFull = TRUE)
  r2 <- .foceiCovAnalytic(fit, sens = "fd2",    covFull = TRUE)
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
  r <- .foceiCov(fit, covFull = TRUE)
  expect_false(is.null(r))
  expect_match(r$method, "^analytic")                # block Omega via the E-basis derivatives
  expect_true(any(grepl("^cov\\.", r$params)))       # an off-diagonal Omega covariance SE
  expect_true(all(is.finite(r$se)) && all(r$se > 0))
  expect_true(all(eigen(r$cov, symmetric = TRUE, only.values = TRUE)$values > 0))
  r2 <- .foceiCovAnalytic(fit, sens = "fd2", covFull = TRUE)
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
  r <- .foceiCovAnalytic(fit, covFull = TRUE)
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
  r <- .foceiCovAnalytic(fa, covFull = TRUE)
  expect_false(is.null(r))
  expect_true("wt_eff" %in% r$params)                # covariate coefficient present
  expect_true(all(is.finite(r$se)) && all(r$se > 0))
  # the fd2 tier now shares the direction-set model, so it handles the covariate
  # direction too (previously it bowed out) and reproduces exact3
  rf <- .foceiCovAnalytic(fa, sens = "fd2", covFull = TRUE)
  expect_false(is.null(rf))
  expect_setequal(rf$params, r$params)
  expect_equal(unname(rf$se[r$params]), unname(r$se), tolerance = 1e-3)
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
  r <- .foceiCovAnalytic(fa, covFull = TRUE)
  expect_false(is.null(r))
  expect_true("tv" %in% r$params)                    # eta-less theta present (would bow out under scaling)
  expect_true(all(is.finite(r$se)) && all(r$se > 0))
  # fd2 handles the eta-less direction too and reproduces exact3
  rf <- .foceiCovAnalytic(fa, sens = "fd2", covFull = TRUE)
  expect_false(is.null(rf))
  expect_setequal(rf$params, r$params)
  expect_equal(unname(rf$se[r$params]), unname(r$se), tolerance = 1e-3)
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
  # the parameter table is wired from the analytic cov even though the FD cov step
  # was skipped: covMethod is labelled and the structural-theta SEs are populated
  expect_equal(fa$covMethod, "analytic-fd2")         # fd2 is the default tier
  .se <- fa$parFixedDf[c("tka", "tcl", "tv"), "SE"]
  expect_true(all(is.finite(.se)) && all(.se > 0))
  ffd <- suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "r", covType = "fd",
                                              covFull = TRUE, sigdig = 6)))
  seA <- sqrt(abs(diag(fa$cov))); seF <- sqrt(abs(diag(ffd$cov)))
  expect_equal(unname(seA[names(seF)]), unname(seF), tolerance = 5e-3)
  # the analytic parFixed SEs equal the FD parFixed SEs (both wired the same way)
  expect_equal(unname(fa$parFixedDf[c("tka","tcl","tv"), "SE"]),
               unname(ffd$parFixedDf[c("tka","tcl","tv"), "SE"]), tolerance = 5e-3)
})

test_that("a non-mu-referenced eta is handled by the uniform direction engine", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  # eta.ka sits on a FIXED population ka (no tka) -> muRefDataFrame does not pair it
  # with a theta.  The engine keeps its ETA_i_ direction and names its Omega by the
  # eta (om.eta.ka); no structural theta reuses that direction.
  nmr <- function() {
    ini({ tcl <- log(2.7); tv <- log(31.5); eta.ka ~ 0.3; eta.cl ~ 0.3; add.sd <- 0.7 })
    model({ ka <- 1.5 * exp(eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) })
  }
  fit <- suppressMessages(nlmixr(nmr, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "", sigdig = 6)))
  r <- .foceiCovAnalytic(fit, covFull = TRUE)
  expect_false(is.null(r))
  expect_true("om.eta.ka" %in% r$params)             # non-mu-ref eta's Omega named by the eta
  expect_true(all(is.finite(r$se)) && all(r$se > 0))
  # fd2 handles it too, and both match the finite-difference full covariance
  rf <- .foceiCovAnalytic(fit, sens = "fd2", covFull = TRUE)
  expect_false(is.null(rf))
  expect_equal(unname(rf$se[r$params]), unname(r$se), tolerance = 2e-3)
  ffd <- suppressMessages(nlmixr(nmr, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "r", covFull = TRUE, sigdig = 6)))
  seF <- sqrt(abs(diag(ffd$cov)))
  expect_equal(unname(r$se[names(seF)]), unname(seF), tolerance = 1e-2)
})

test_that("covType='analytic' falls back to the finite-difference cov out of scope", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  # the lnorm error model is out of analytic scope, so the engine announces it and
  # the finite-difference covariance runs as the fallback
  oos <- function() {
    ini({
      tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; ln.sd <- 0.3
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ lnorm(ln.sd)
    })
  }
  fit <- suppressMessages(suppressWarnings(
    nlmixr(oos, nlmixr2data::theo_sd, "focei",
           foceiControl(print = 0L, covMethod = "r", covType = "analytic", sigdig = 6))))
  expect_false(grepl("analytic", fit$covMethod))       # FD ran as the fallback
  expect_true(is.matrix(fit$cov))
})

test_that("the combined1 error model is handled and matches the finite-difference cov", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  cmb <- function() {
    ini({
      tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.3; prop.sd <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) + prop(prop.sd)
    })
  }
  fa <- suppressMessages(nlmixr(cmb, nlmixr2data::theo_sd, "focei",
                                foceiControl(print = 0L, covMethod = "", addProp = "combined1", sigdig = 6)))
  r <- .foceiCovAnalytic(fa, covFull = TRUE)             # combined1: R = (sa + sp f)^2
  expect_false(is.null(r))
  expect_true(all(is.finite(r$se)) && all(r$se > 0))
  ffd <- suppressMessages(nlmixr(cmb, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "r", covFull = TRUE,
                                              addProp = "combined1", sigdig = 6)))
  seF <- sqrt(abs(diag(ffd$cov)))
  expect_equal(unname(r$se[names(seF)]), unname(seF), tolerance = 1e-2)
})

test_that("analytic covariance is FOCEI only; FOCE falls back to finite differences", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  # The analytic engine is scoped to FOCEI.  FOCE (no interaction) has a non-smooth
  # objective (its inner EBE solve is fed an inconsistent value/gradient pair), so it
  # is left to its own finite-difference cov path (tracked separately in #694).
  prp <- function() {
    ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5); eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.2; prop.sd <- 0.1 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) + prop(prop.sd) })
  }
  fpi <- suppressMessages(nlmixr(prp, nlmixr2data::theo_sd, "focei", foceiControl(print = 0L, covMethod = "", sigdig = 6)))
  fpe <- suppressMessages(nlmixr(prp, nlmixr2data::theo_sd, "foce",  foceiControl(print = 0L, covMethod = "", sigdig = 6)))
  expect_true(nlmixr2est:::.foceiAnalyticInScope(fpi$finalUi))    # FOCEI in scope
  expect_false(nlmixr2est:::.foceiAnalyticInScope(fpe$finalUi))   # FOCE out of scope -> FD fallback
  rpi <- .foceiCovAnalytic(fpi, covFull = TRUE)
  expect_false(is.null(rpi))                                      # FOCEI: analytic cov
  expect_null(.foceiCovAnalytic(fpe, covFull = TRUE))             # FOCE: bows out
})
