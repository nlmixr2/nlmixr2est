# Analytic FOCEI covariance: the standalone engine (foceiCov / foceiCovAnalytic)
# and the covMethod="analytic" seam (analytic R-matrix while the optimizer is live,
# with a finite-difference fallback out of scope).

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
  skip_on_ci()
  fit <- suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "")))
  r <- foceiCovAnalytic(fit)
  expect_false(is.null(r))
  expect_identical(r$method, "analytic")
  # every population parameter present: 3 theta + 1 sigma + 3 Omega variances
  expect_setequal(r$params, c("tka", "tcl", "tv", "add.sd", "om.eta.ka", "om.eta.cl", "om.eta.v"))
  expect_true(all(is.finite(r$se)) && all(r$se > 0))
  # structural-theta SEs match NONMEM $COV MATRIX=R (0.19180/0.08352/0.04661)
  expect_equal(unname(r$se[c("tka", "tcl", "tv")]),
               c(0.18868, 0.08351, 0.04617), tolerance = 0.03)
  # R matrix symmetric positive-definite
  expect_lt(max(abs(r$R - t(r$R))), 1e-6)
  expect_true(all(eigen(r$R, symmetric = TRUE, only.values = TRUE)$values > 0))
})

test_that("block Omega is handled analytically, with off-diagonal covariance SEs", {
  skip_on_cran()
  skip_on_ci()
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
  # theo_sd alone (12 subjects) does not identify the off-diagonal covariance:
  # its full observed information is genuinely indefinite (the exact analytic R and
  # a finite-difference R agree on a small negative eigenvalue), so cov(eta.cl,eta.v)
  # has a negative variance and a NaN SE under BOTH methods.  Replicate the data so
  # the block is estimable and the full covariance is positive-definite -- then the
  # analytic engine yields a real off-diagonal covariance SE.
  d0 <- nlmixr2data::theo_sd
  dat <- do.call(rbind, lapply(1:4, function(k) { .x <- d0; .x$ID <- .x$ID + (k - 1) * 100; .x }))
  fit <- suppressMessages(nlmixr(blk, dat, "focei",
                                 foceiControl(print = 0L, covMethod = "")))
  r <- foceiCovAnalytic(fit)
  expect_false(is.null(r))
  expect_identical(r$method, "analytic")             # block Omega via the E-basis derivatives
  expect_true(any(grepl("^cov\\.", r$params)))       # an off-diagonal Omega covariance SE
  expect_true(all(is.finite(r$se)) && all(r$se > 0))
  expect_true(all(eigen(r$cov, symmetric = TRUE, only.values = TRUE)$values > 0))
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
  r <- foceiCovAnalytic(fit)
  expect_false(is.null(r))
  expect_setequal(r$params, c("tcl", "add.sd", "om.eta.cl"))
  expect_true(all(is.finite(r$se)) && all(r$se > 0))
})

test_that("covMethod='analytic' installs the full analytic covariance on the fit", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  fit <- suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "analytic", covFull = TRUE)))
  # covFull=TRUE swaps in the full theta+sigma+Omega cov (7x7), not the theta-only FD cov
  expect_true(is.matrix(fit$cov))
  expect_setequal(rownames(fit$cov),
                  c("tka", "tcl", "tv", "add.sd", "om.eta.ka", "om.eta.cl", "om.eta.v"))
  .se <- sqrt(diag(fit$cov))
  expect_equal(unname(.se[c("tka", "tcl", "tv")]),
               c(0.18868, 0.08351, 0.04617), tolerance = 0.03)
  expect_true(all(is.finite(.se)) && all(.se > 0))
})

test_that("finite-difference covMethod='r,s' covFull=TRUE installs the full FD covariance matching analytic", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  # the finite-difference covariance over the SAME full parameter set as the analytic
  # engine: structural + residual thetas plus the Omega variance-covariance elements
  # (Gill-style adaptive step; Omega perturbed on the variance scale, no Jacobian).
  fa <- suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
                                foceiControl(print = 0L, covMethod = "analytic", covFull = TRUE)))
  ff <- suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
                                foceiControl(print = 0L, covMethod = "r,s", covFull = TRUE)))
  # full theta+sigma+Omega cov, same parameter set as analytic
  expect_setequal(rownames(ff$cov),
                  c("tka", "tcl", "tv", "add.sd", "om.eta.ka", "om.eta.cl", "om.eta.v"))
  .seF <- sqrt(diag(ff$cov)); .seA <- sqrt(diag(fa$cov))[rownames(ff$cov)]
  expect_true(all(is.finite(.seF)) && all(.seF > 0))
  # FD full SEs match the analytic full SEs (finite-difference tolerance)
  expect_equal(unname(.seF), unname(.seA), tolerance = 0.05)
})

test_that("finite-difference covMethod='r,s' covFull=FALSE keeps the finite-difference theta covariance", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  fit <- suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
                                 foceiControl(print = 0L, covMethod = "r,s", covFull = FALSE)))
  # FD covFull=FALSE cov is theta-only (residual/Omega are skipCov'd): no Omega/residual
  # rows, and it is a valid, finite, positive covariance -- not the full analytic matrix
  expect_true(is.matrix(fit$cov))
  expect_false(any(grepl("^om\\.", rownames(fit$cov))))
  expect_true(all(is.finite(diag(fit$cov))) && all(diag(fit$cov) > 0))
})

test_that("covMethod='analytic' falls back to the finite-difference cov out of scope", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  # the laplace censored determinant is out of analytic-covariance scope -> analytic bails and
  # the live finite-difference Hessian is used (a valid theta cov).  (censOption="gauss"
  # censored IS in scope now -- see the "covers censored M2/M3/M4" test below.)
  cm <- function() {
    ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
          eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) })
  }
  d <- nlmixr2data::theo_sd
  d$CENS <- ifelse(d$DV < 2 & d$EVID == 0, 1L, 0L); d$DV[d$CENS == 1] <- 2
  # out of scope -> foceiCalcR warns (visibly) and uses the finite-difference cov
  fit <- suppressWarnings(suppressMessages(nlmixr(cm, d, "focei",
                                                  foceiControl(print = 0L, covMethod = "analytic", censOption = "laplace"))))
  expect_true(is.matrix(fit$cov))
  expect_false(any(grepl("^om\\.", rownames(fit$cov))))  # theta-only FD cov, not the full analytic
})

test_that("covMethod='analytic' covers censored M2/M3/M4 for FOCEI and FOCE (gauss); laplace uses FD", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  cm <- function() {
    ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
          eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) })
  }
  base <- nlmixr2data::theo_sd
  dM3 <- base; dM3$CENS <- ifelse(dM3$DV < 2 & dM3$EVID == 0, 1L, 0L); dM3$DV[dM3$CENS == 1] <- 2
  dM2 <- base; dM2$CENS <- 0L; dM2$LIMIT <- 0
  dM4 <- dM3; dM4$LIMIT <- 0
  # FOCEI + gauss: full analytic cov (theta+sigma+Omega), close to the FD Hessian cov
  fitA <- suppressWarnings(suppressMessages(nlmixr(cm, dM3, "focei",
                                                   foceiControl(print = 0L, covMethod = "analytic", sigdig = 6))))
  expect_identical(fitA$covMethod, "analytic")
  expect_true(any(grepl("^om\\.", rownames(fitA$cov))))
  expect_true(all(is.finite(sqrt(diag(fitA$cov)))))
  fitR <- suppressWarnings(suppressMessages(nlmixr(cm, dM3, "focei",
                                                   foceiControl(print = 0L, covMethod = "r", sigdig = 6))))
  cp <- intersect(rownames(fitA$cov), rownames(fitR$cov))
  expect_lt(max(abs(sqrt(diag(fitA$cov))[cp] - sqrt(diag(fitR$cov))[cp]) /
                  (sqrt(diag(fitR$cov))[cp] + 1e-8)), 0.05)
  for (dd in list(dM2, dM4)) {
    f <- suppressWarnings(suppressMessages(nlmixr(cm, dd, "focei",
                                                  foceiControl(print = 0L, covMethod = "analytic"))))
    expect_identical(f$covMethod, "analytic")
    expect_true(any(grepl("^om\\.", rownames(f$cov))))
  }
  # FOCE (gauss) censored is in scope too: full analytic cov, theta/sigma SEs close to FD
  fF <- suppressWarnings(suppressMessages(nlmixr(cm, dM3, "focei",
                                                 foceiControl(print = 0L, covMethod = "analytic", interaction = FALSE, sigdig = 6))))
  expect_identical(fF$covMethod, "analytic")
  expect_true(any(grepl("^om\\.", rownames(fF$cov))))
  fFr <- suppressWarnings(suppressMessages(nlmixr(cm, dM3, "focei",
                                                  foceiControl(print = 0L, covMethod = "r", interaction = FALSE, sigdig = 6))))
  cpf <- intersect(c("tka", "tcl", "tv", "add.sd"), intersect(rownames(fF$cov), rownames(fFr$cov)))
  expect_lt(max(abs(sqrt(diag(fF$cov))[cpf] - sqrt(diag(fFr$cov))[cpf]) /
                  (sqrt(diag(fFr$cov))[cpf] + 1e-8)), 0.03)
  # foce+ (live conditional R) censored is in scope too
  fFp <- suppressWarnings(suppressMessages(nlmixr(cm, dM3, "focei",
                                                  foceiControl(print = 0L, covMethod = "analytic",
                                                               interaction = FALSE, foceType = "foce+"))))
  expect_identical(fFp$covMethod, "analytic")
  expect_true(any(grepl("^om\\.", rownames(fFp$cov))))
  # only the laplace censored determinant stays on the FD cov (theta-only)
  fL <- suppressWarnings(suppressMessages(nlmixr(cm, dM3, "focei",
                                                 foceiControl(print = 0L, covMethod = "analytic", censOption = "laplace"))))
  expect_false(any(grepl("^om\\.", rownames(fL$cov))))
})

test_that("covMethod='analytic' with pure proportional error near a zero prediction falls back to FD", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  # pure proportional error IS in analytic scope, but its variance sp^2 f^2 vanishes as
  # f -> 0.  theo_sd is oral, so the predicted concentration is ~0 at the pre-dose time:
  # the near-zero-prediction guard must catch it and give a valid FD cov, never crash.
  pm <- function() {
    ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
          eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; prop.sd <- 0.2 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ prop(prop.sd) })
  }
  fit <- suppressWarnings(suppressMessages(nlmixr(pm, nlmixr2data::theo_sd, "focei",
                                                  foceiControl(print = 0L, covMethod = "analytic"))))
  expect_true(is.matrix(fit$cov))
  expect_false(any(grepl("^om\\.", rownames(fit$cov))))
})

test_that("foce+ (live-R) additive analytic R equals the FOCEI analytic R at the same theta", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  # additive error: R = sa^2 is constant, so live vs frozen R and the FOCEI interaction
  # term all coincide -- the foce+ analytic R must reproduce FOCEI's.  maxOuterIterations=0
  # evaluates both at the identical initial theta, so the comparison is tight (the only
  # slack is the inner EBE tolerance).
  fitP <- suppressWarnings(suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
            foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L,
                         interaction = FALSE, foce = "foce+"))))
  fitI <- suppressWarnings(suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
            foceiControl(print = 0L, covMethod = "", maxOuterIterations = 0L))))
  rP <- foceiCovAnalytic(fitP); rI <- foceiCovAnalytic(fitI)
  expect_false(is.null(rP)); expect_identical(rP$method, "analytic")
  expect_false(is.null(rI))
  expect_lt(max(abs(rP$R - rI$R) / (abs(rI$R) + 1e-8)), 1e-3)
})

# Wang 2007 monoexponential IV bolus: predictions 10*exp(-ke*t) are bounded away from
# zero at every observation, so pure proportional error is genuinely in analytic scope.
.cov_wang_prop <- function() {
  ini({ tke <- log(0.5); eta.ke ~ 0.04; prop.sd <- sqrt(0.1) })
  model({ ke <- exp(tke + eta.ke); d/dt(ipre) <- -ke * ipre; ipre ~ prop(prop.sd) })
}
.cov_wang_data <- function() {
  d <- nlmixr2data::Wang2007; d$DV <- d$Y
  dose <- d[d$Time == 0, ]; dose$EVID <- 101; dose$AMT <- 10
  dat <- rbind(dose, data.frame(d, EVID = 0, AMT = 0))
  dat[order(dat$ID, -dat$EVID, dat$Time), ]
}

test_that("covMethod='analytic' handles pure proportional error away from zero (FOCEI and FOCE)", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  dat <- .cov_wang_data()
  for (est in c("focei", "foce")) {
    fit <- suppressMessages(nlmixr(.cov_wang_prop, dat, est,
                                   foceiControl(print = 0L, covMethod = "analytic", covFull = TRUE)))
    r <- foceiCovAnalytic(fit)
    expect_false(is.null(r))                                   # in scope, not an FD fallback
    expect_identical(r$method, "analytic")
    expect_setequal(r$params, c("tke", "prop.sd", "om.eta.ke"))
    expect_true(all(is.finite(r$se)) && all(r$se > 0))
    # the full analytic cov is installed on the fit (om. row present)
    expect_true(any(grepl("^om\\.", rownames(fit$cov))))
    # analytic SEs match a Richardson finite-difference of the objective (validated to
    # ~1e-4 vs numDeriv); the reference values are the converged plateau / NONMEM MATRIX=R
    .ref <- if (est == "focei") c(tke = 0.09234, prop.sd = 0.007446, om.eta.ke = 0.03684)
            else                c(tke = 0.09065, prop.sd = 0.007624, om.eta.ke = 0.03624)
    expect_equal(unname(r$se[names(.ref)]), unname(.ref), tolerance = 0.01)
  }
})

test_that("estimated boxCox lambda: analytic cov (FOCEI/FOCE/foce+) matches the s estimator", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("nlmixr2data")
  mBox <- function() {
    ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45; eta.ka ~ 0.5; eta.cl ~ 0.08; eta.v ~ 0.05
          add.sd <- 0.7; lambda <- c(-2, 0.9, 3) })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
            d/dt(depot) <- -ka * depot; d/dt(center) <- ka * depot - cl / v * center
            cp <- center / v; cp ~ add(add.sd) + boxCox(lambda) })
  }
  d <- nlmixr2data::theo_sd
  chk <- function(est, ctlExtra = list()) {
    ctlA <- do.call(foceiControl, c(list(print = 0L, covMethod = "analytic", covFull = TRUE, fast = TRUE), ctlExtra))
    ctlS <- do.call(foceiControl, c(list(print = 0L, covMethod = "s", covFull = TRUE, fast = TRUE), ctlExtra))
    fitA <- suppressMessages(nlmixr2(mBox, d, est, ctlA))
    fitS <- suppressMessages(nlmixr2(mBox, d, est, ctlS))
    expect_identical(fitA$covMethod, "analytic")       # analytic ran (not an FD fallback)
    seA <- sqrt(diag(fitA$cov)); seS <- sqrt(diag(fitS$cov))
    nm <- c("tka", "tcl", "tv", "add.sd", "lambda")     # theta/sigma/lambda block (DV-affected)
    expect_true(all(is.finite(seA[nm])) && all(seA[nm] > 0))
    expect_equal(unname(seA[nm]), unname(seS[nm]), tolerance = 0.05)
  }
  chk("focei")
  chk("foce")
  chk("foce", list(foce = "foce+"))                    # focep: residual at the posthoc eta
})

test_that("covMethod='analytic' emits an informative message when it falls back to FD", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  # the laplace censored determinant is out of scope (only the gauss censored cov is ported);
  # with covMethod="analytic" the fallback to the FD cov is announced (message, not warning)
  # so the user knows why they did not get analytic.  (lnorm, fixed/estimated-lambda transforms
  # and the default gauss censored cov ARE in scope.)
  cm <- function() {
    ini({ tcl <- log(2.7); eta.cl ~ 0.1; add.sd <- 0.7 })
    model({ ka <- 1.5; cl <- exp(tcl + eta.cl); v <- 31.5
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) })
  }
  d <- nlmixr2data::theo_sd
  d$CENS <- ifelse(d$DV < 2 & d$EVID == 0, 1L, 0L); d$DV[d$CENS == 1] <- 2
  expect_message(
    suppressWarnings(nlmixr(cm, d, "focei",
                            foceiControl(print = 0L, covMethod = "analytic", censOption = "laplace"))),
    "covType=\"analytic\".*finite-difference")
})

test_that("covMethod='analytic' handles a non-mu-referenced covariate coefficient", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  # a non-mu-referenced covariate coefficient (wt_cl) gets its own THETA-direction
  # sensitivity: the analytic full theta+sigma+Omega cov is installed, and the
  # covariate SE is finite/positive and matches a tight-tolerance finite-difference.
  cvm <- function() {
    ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5); wt_cl <- 0.75
          eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + wt_cl * log(WT / 70) + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) })
  }
  dat <- nlmixr2data::theo_sd
  fit <- suppressWarnings(suppressMessages(nlmixr(cvm, dat, "focei",
                                                  foceiControl(print = 0L, covMethod = "analytic", covFull = TRUE))))
  expect_true(is.matrix(fit$cov))
  # covFull=TRUE installs the full analytic cov: covariate theta + Omega rows are present
  expect_true("wt_cl" %in% rownames(fit$cov))
  expect_true(any(grepl("^om\\.", rownames(fit$cov))))
  .se <- sqrt(diag(fit$cov))
  # covariate SE is finite and positive (not the mu-ref bail)
  expect_true(is.finite(.se[["wt_cl"]]) && .se[["wt_cl"]] > 0)
  # mu-ref thetas still match NONMEM $COV MATRIX=R (0.18868/0.08351/0.04617)
  expect_equal(unname(.se[c("tka", "tcl", "tv")]),
               c(0.18868, 0.08351, 0.04617), tolerance = 0.05)
  # covariate SE matches a tight-tolerance finite-difference R covariance (~0.607)
  fit_fd <- suppressWarnings(suppressMessages(nlmixr(cvm, dat, "focei",
              foceiControl(print = 0L, covMethod = "r", sigdig = 7))))
  .sefd <- sqrt(diag(fit_fd$cov))
  expect_equal(unname(.se[["wt_cl"]]), unname(.sefd[["wt_cl"]]), tolerance = 0.05)
})

test_that("covMethod='analytic' handles a non-mu-referenced eta (orphan Omega variance)", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  # cl <- exp(tcl)*exp(eta.cl) is NOT recognized as mu-referenced by nlmixr, so eta.cl
  # is an orphan eta: it keeps its own ETA-direction and an Omega variance named by the
  # eta (om.eta.cl).  The full analytic cov must install and the theta SEs must match a
  # tight-tolerance finite-difference R covariance.
  nonmu <- function() {
    ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
          eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka)
      cl <- exp(tcl) * exp(eta.cl)                 # non-mu-referenced clearance
      v  <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) })
  }
  dat <- nlmixr2data::theo_sd
  fit <- suppressWarnings(suppressMessages(nlmixr(nonmu, dat, "focei",
                                                  foceiControl(print = 0L, covMethod = "analytic", covFull = TRUE))))
  expect_true(is.matrix(fit$cov))
  # the orphan eta's variance is named by the eta, not a theta
  expect_true("om.eta.cl" %in% rownames(fit$cov))
  .se <- sqrt(diag(fit$cov))
  expect_true(is.finite(.se[["om.eta.cl"]]) && .se[["om.eta.cl"]] > 0)
  # theta SEs match a tight-tolerance finite-difference R covariance
  fit_fd <- suppressWarnings(suppressMessages(nlmixr(nonmu, dat, "focei",
              foceiControl(print = 0L, covMethod = "r", sigdig = 7))))
  .sefd <- sqrt(diag(fit_fd$cov))
  expect_equal(unname(.se[c("tka", "tcl", "tv")]),
               unname(.sefd[c("tka", "tcl", "tv")]), tolerance = 0.05)
})

test_that("covFull controls fit$cov shape without changing the theta SEs", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  # covFull=TRUE installs the full theta+sigma+Omega cov; covFull=FALSE installs the
  # theta block (structural + residual, i.e. the non-skipped thetas), no Omega.  The
  # theta SEs are identical either way.
  fitT <- suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
                                  foceiControl(print = 0L, covMethod = "analytic", covFull = TRUE)))
  fitF <- suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
                                  foceiControl(print = 0L, covMethod = "analytic", covFull = FALSE)))
  fitD <- suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
                                  foceiControl(print = 0L, covMethod = "analytic")))
  .th <- c("tka", "tcl", "tv", "add.sd")   # non-skipped thetas: structural + residual
  # covFull=TRUE adds the Omega block; covFull=FALSE is the theta block (no Omega)
  expect_true(any(grepl("^om\\.", rownames(fitT$cov))))
  expect_false(any(grepl("^om\\.", rownames(fitF$cov))))
  expect_setequal(rownames(fitF$cov), .th)
  expect_true(any(grepl("^om\\.", rownames(fitD$cov))))   # default is now covFull=TRUE (the full cov)
  # identical theta SEs, and covFull=FALSE is exactly the theta submatrix of the full cov
  expect_equal(sqrt(diag(fitF$cov))[.th], sqrt(diag(fitT$cov))[.th])
  expect_equal(unname(fitF$cov[.th, .th]), unname(fitT$cov[.th, .th]))
})

test_that("covMethod='analytic' restores the fit solve so tables stay intact", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  d <- nlmixr2data::theo_sd
  # the augmented sensitivity solves replace the global solve; it must be restored
  # or foceiFinalizeTables reads the last subject's solve (truncated per-obs tables)
  fa <- suppressMessages(nlmixr(.cov_one_cmt, d, "focei",
                                foceiControl(print = 0L, covMethod = "analytic")))
  ff <- suppressMessages(nlmixr(.cov_one_cmt, d, "focei",
                                foceiControl(print = 0L, covMethod = "r,s")))
  expect_equal(nrow(fa), nrow(ff))
})

test_that(".omegaBlocks uses the declared block, not converged values", {
  # a declared 2x2 block (etas 2-3) whose off-diagonal converges near zero must stay
  # one block, so its covariance parameter is not silently dropped from the cov
  Om <- diag(3); Om[2, 3] <- Om[3, 2] <- 1e-12
  idf <- data.frame(neta1 = c(1, 2, 3, 2), neta2 = c(1, 2, 3, 3))
  blk <- nlmixr2est:::.omegaBlocks(Om, idf) # nolint: undesirable_operator_linter.
  expect_length(blk, 2L)
  expect_true(any(vapply(blk, function(b) all(c(2L, 3L) %in% b), logical(1))))
})

test_that("covMethod='analytic' joins subjects by ID code, not factor label (non-1..N IDs)", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  # subjects keyed by the etTrans integer code, not the original ID label: relabeling
  # theo_sd's IDs to 101..112 (or permuting them) must not change the covariance.
  d1 <- nlmixr2data::theo_sd                 # IDs 1..12
  d2 <- d1; d2$ID <- d2$ID + 100L            # IDs 101..112 (non-1..N)
  set.seed(1); pm <- sample(1:12); d3 <- d1; d3$ID <- pm[d1$ID]   # a permutation of 1..N
  ctl <- foceiControl(print = 0L, covMethod = "analytic", covFull = TRUE)
  f1 <- suppressMessages(nlmixr(.cov_one_cmt, d1, "focei", ctl))
  f2 <- suppressMessages(nlmixr(.cov_one_cmt, d2, "focei", ctl))
  f3 <- suppressMessages(nlmixr(.cov_one_cmt, d3, "focei", ctl))
  expect_true(any(grepl("^om\\.", rownames(f2$cov))))          # analytic ran (not silent FD)
  expect_setequal(rownames(f1$cov), rownames(f2$cov))
  # SEs identical to the 1..N fit (a wrong join would silently pair the wrong events)
  expect_equal(sqrt(diag(f2$cov)), sqrt(diag(f1$cov)))
  expect_equal(sqrt(diag(f3$cov)), sqrt(diag(f1$cov)))
})

test_that("covMethod='analytic' falls back to FD when a theta is shared by two etas", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  # tcl mu-references BOTH eta.cl and eta.v: the analytic direction map cannot send one
  # theta down two eta routes, so it must bow out to the (correct) finite-difference cov.
  twoEta <- function() {
    ini({ tka <- log(1.5); tcl <- log(2.7)
          eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tcl + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v; cp ~ add(add.sd) })
  }
  fit <- suppressWarnings(suppressMessages(nlmixr(twoEta, nlmixr2data::theo_sd, "focei",
                          foceiControl(print = 0L, covMethod = "analytic", covFull = TRUE))))
  expect_true(is.matrix(fit$cov))
  expect_false(any(grepl("^om\\.", rownames(fit$cov))))        # theta-only FD cov, not the full analytic
})

test_that("covMethod='analytic' falls back to FD under a bounded-parameter transform", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  skip_if_not_installed("minqa")
  # a finitely-bounded theta + an outer optimizer without native bounds (newuoa) rewrites
  # the model to an internal scale; the Jacobian hook corrects env$cov.  Analytic must bow
  # out (else it overwrites the Jacobian-corrected cov with the internal-scale one).
  bnd <- function() {
    ini({ tka <- log(1.5); tcl <- log(2.7); tv <- c(2, log(31.5), 5)
          eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v; cp ~ add(add.sd) })
  }
  fit <- suppressWarnings(suppressMessages(nlmixr(bnd, nlmixr2data::theo_sd, "focei",
                          foceiControl(print = 0L, covMethod = "analytic", covFull = TRUE, outerOpt = "newuoa"))))
  expect_true(is.matrix(fit$cov))
  expect_false(any(grepl("^om\\.", rownames(fit$cov))))        # FD (Jacobian-correct), not analytic
})

test_that("covMethod='analytic' handles SD-scale IOV and falls back for other iovXform", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  # simulate a small data set with genuine between-occasion variability in CL so IOV is
  # identified (the fixture theo data has none, driving iov.cl onto its boundary).
  set.seed(42)
  .sim <- rxode2::rxode2({
    ka <- exp(0.5 + eta.ka); cl <- exp(1.0 + eta.cl + iov.cl); v <- exp(3.4 + eta.v)
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - cl / v * center
    cp <- center / v
  })
  .ev <- rxode2::et(seq(0.25, 20, by = 2)); .ev <- rxode2::et(.ev, seq(48.25, 68, by = 2))
  .ev <- rxode2::et(.ev, amt = 320, time = 0, cmt = "depot")
  .ev <- rxode2::et(.ev, amt = 320, time = 48, cmt = "depot")
  .ev <- as.data.frame(.ev)
  .nid <- 24
  .d <- do.call(rbind, lapply(seq_len(.nid), function(id) {
    .dd <- .ev
    .dd$iov.cl <- ifelse(.dd$time < 40, rnorm(1, 0, sqrt(0.08)), rnorm(1, 0, sqrt(0.08)))
    .s <- rxode2::rxSolve(.sim, c(eta.ka = rnorm(1, 0, sqrt(0.3)), eta.cl = rnorm(1, 0, sqrt(0.1)),
                                  eta.v = rnorm(1, 0, sqrt(0.1))), .dd, returnType = "data.frame")
    .obs <- .s[!is.na(.s$cp) & .s$cp > 0, ]
    data.frame(ID = id, TIME = .obs$time, DV = .obs$cp * exp(rnorm(nrow(.obs), 0, 0.1)),
               AMT = 0, EVID = 0, occ = ifelse(.obs$time < 40, 1L, 2L))
  }))
  .dose <- data.frame(ID = rep(seq_len(.nid), each = 2), TIME = rep(c(0, 48), .nid), DV = 0,
                      AMT = 320, EVID = 1, occ = rep(c(1L, 2L), .nid))
  dat <- rbind(.dose, .d); dat <- dat[order(dat$ID, dat$TIME, -dat$EVID), ]
  iovm <- function() {
    ini({ tka <- 0.5; tcl <- 1.0; tv <- 3.4
          eta.ka ~ 0.3; eta.cl ~ 0.1; eta.v ~ 0.1
          iov.cl ~ 0.08 | occ
          add.sd <- 0.3 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl + iov.cl); v <- exp(tv + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v; cp ~ add(add.sd) })
  }
  # default sd scale: analytic path installs the full cov with the IOV variance row
  fSD <- suppressWarnings(suppressMessages(nlmixr(iovm, dat, "focei",
              foceiControl(print = 0L, covMethod = "analytic", covFull = TRUE, iovXform = "sd"))))
  expect_true(is.matrix(fSD$cov))
  expect_true(any(grepl("^om\\.", rownames(fSD$cov))))         # analytic ran
  expect_true("iov.cl" %in% rownames(fSD$cov))                 # IOV shared-variance SE present
  # non-sd iovXform uses a different predictor/chain-rule -> must fall back to FD
  fVAR <- suppressWarnings(suppressMessages(nlmixr(iovm, dat, "focei",
              foceiControl(print = 0L, covMethod = "analytic", covFull = TRUE, iovXform = "var"))))
  expect_true(is.matrix(fVAR$cov))
  expect_false(any(grepl("^om\\.", rownames(fVAR$cov))))       # theta-only FD cov
})

test_that("covMethod selects the analytic-vs-FD seam and the reporting formula", {
  # covType was folded into covMethod: "analytic" is the exact observed-information R
  # (reported with the "r" formula, so covMethod=2L) carried to the solver as the internal
  # covType="analytic"; the finite-difference formulas keep covType="fd"; "" skips cov.
  .ca <- foceiControl(covMethod = "analytic")
  expect_identical(.ca$covMethod, 2L)
  expect_identical(.ca$covType, "analytic")
  expect_identical(foceiControl(covMethod = "r,s")$covMethod, 1L)
  expect_identical(foceiControl(covMethod = "r")$covMethod, 2L)
  expect_identical(foceiControl(covMethod = "s")$covMethod, 3L)
  expect_identical(foceiControl(covMethod = "r")$covType, "fd")
  expect_identical(foceiControl(covMethod = "")$covMethod, 0L)  # "" skips the covariance step
  # analytic is the default
  .cd <- foceiControl()
  expect_identical(.cd$covMethod, 2L)
  expect_identical(.cd$covType, "analytic")
})

test_that("covMethod='analytic' covFull=FALSE respects skipCov (matches the FD shape)", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  # skipCov excludes tv from the theta cov: covFull=FALSE must install the same (2-theta)
  # shape as the finite-difference covMethod, not widen back to every structural theta.
  fitA <- suppressWarnings(suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
              foceiControl(print = 0L, covMethod = "analytic", covFull = FALSE, skipCov = c(FALSE, FALSE, TRUE, TRUE)))))
  fitF <- suppressWarnings(suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
              foceiControl(print = 0L, covMethod = "r", covFull = FALSE, skipCov = c(FALSE, FALSE, TRUE, TRUE)))))
  expect_setequal(rownames(fitA$cov), rownames(fitF$cov))
  expect_false("tv" %in% rownames(fitA$cov))                   # skipCov'd theta excluded, not widened
})

.cov_combined <- function() {
  ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.2; prop.sd <- 0.1 })
  model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - cl / v * center
    cp <- center / v; cp ~ add(add.sd) + prop(prop.sd) })
}

test_that("FOCE (interaction=FALSE) additive analytic cov equals the FOCEI analytic cov", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  # additive error: R does not depend on eta, so the FOCEI interaction term (dR/deta) is
  # identically 0 and the FOCE (interaction=0) analytic covariance coincides with FOCEI.
  fF <- suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
          foceiControl(print = 0L, covMethod = "analytic", covFull = TRUE, interaction = FALSE)))
  fI <- suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focei",
          foceiControl(print = 0L, covMethod = "analytic", covFull = TRUE)))
  expect_true(any(grepl("^om\\.", rownames(fF$cov))))          # analytic ran (not silent FD)
  .th <- c("tka", "tcl", "tv")
  seF <- sqrt(diag(fF$cov))[.th]; seI <- sqrt(diag(fI$cov))[.th]
  expect_true(all(is.finite(seF)) && all(seF > 0))
  expect_equal(unname(seF), unname(seI), tolerance = 0.01)
})

test_that("FOCE (interaction=FALSE) combined analytic cov matches the corrected-FOCE gold FD", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  # CORRECTED FOCE: the estimator freezes the residual variance R0 at the eta=0
  # POPULATION prediction (getPopR); the analytic path builds q0=-(y-f)/R0, p=1/R0 and
  # their theta-chain from an eta=0 augmented solve.  Validation is against a FULLY
  # INDEPENDENT finite-difference Hessian of the corrected-FOCE objective
  # (Phi(R0) + 0.5 log|H~_FOCE|, R0 at the eta=0 prediction, EBEs re-solved to
  # S_FOCE = sum(-(y-f)/R0 . a) + Omega^-1 eta = 0).  The FOCE objective differs from
  # FOCEI, so the covariances differ.  NB: for this model+data the corrected FOCE
  # observed information is indefinite (one negative eigenvalue), so a couple of SEs
  # are NaN in BOTH the analytic and the gold FD -- the correctness criterion is that
  # the analytic R MATRIX reproduces the gold-FD Hessian, not that it is invertible.
  theo <- nlmixr2data::theo_sd
  fitF <- suppressMessages(nlmixr(.cov_combined, theo, "focei",
            foceiControl(print = 0L, covMethod = "", interaction = FALSE, sigdig = 6)))
  fitI <- suppressMessages(nlmixr(.cov_combined, theo, "focei",
            foceiControl(print = 0L, covMethod = "", sigdig = 6)))
  rF <- foceiCovAnalytic(fitF); rI <- foceiCovAnalytic(fitI)
  expect_false(is.null(rF)); expect_identical(rF$method, "analytic")
  # FOCE combined != FOCEI combined (the interaction term is non-zero for prop error)
  expect_false(isTRUE(all.equal(unname(rF$R), unname(rI$R), tolerance = 1e-2)))

  # gold standard: central-FD Hessian of the CORRECTED FOCE objective (R0 = eta=0
  # population variance), EBEs re-solved to S_FOCE=0 at each perturbed parameter vector.
  ui <- fitF$finalUi; neta <- 3L; etav <- paste0("ETA_", 1:neta, "_")
  am <- .foceiAnalyticAugModelDirs(ui, etav)
  thNames <- names(fitF$theta)
  thBase <- setNames(as.numeric(fitF$theta[thNames]), paste0("THETA_", seq_along(thNames), "_"))
  iTh <- match(c("tka", "tcl", "tv", "add.sd", "prop.sd"), thNames)
  Om0 <- fitF$omega
  byId <- split(fitF$dataSav, as.character(fitF$dataSav$ID))
  ids <- fitF$eta$ID; idCode <- as.integer(ids)
  eta0m <- as.matrix(fitF$eta[, c("eta.ka", "eta.cl", "eta.v")])
  subj <- lapply(seq_along(ids), function(i) {
    s <- byId[[as.character(idCode[i])]]; obs <- s[s$EVID == 0, , drop = FALSE]
    list(s = s, times = obs$TIME, y = obs$DV, eta0 = eta0m[i, ]) })
  .fa <- function(th, eta, s, times) .foceiAnalyticSolveFA(am, c(th, setNames(eta, etav)), s, times, tol = 1e-12)
  objFOCE <- function(psi) {
    th <- thBase; th[iTh] <- psi[1:5]; sa <- psi[4]; sp <- psi[5]
    Om <- Om0; diag(Om) <- psi[6:8]; Oi <- solve(Om); ldOm <- log(det(Om)); tot <- 0
    for (sj in subj) {
      y <- sj$y; s <- sj$s; times <- sj$times
      E0 <- .fa(th, rep(0, neta), s, times); if (is.null(E0)) return(NA_real_)
      R0 <- sa^2 + sp^2 * E0$f^2                       # eta=0 population variance (fixed in eta)
      eta <- sj$eta0
      for (it in 1:100) {
        E <- .fa(th, eta, s, times); if (is.null(E)) return(NA_real_)
        q0 <- -(y - E$f) / R0
        S <- as.numeric(Oi %*% eta); for (l in 1:neta) S[l] <- S[l] + sum(q0 * E$a[, l])
        if (max(abs(S)) < 1e-12) break
        Hf <- Oi; for (l in 1:neta) for (m in 1:neta) Hf[l, m] <- Hf[l, m] + sum((1/R0) * E$a[, l] * E$a[, m] + q0 * E$A[, l, m])
        eta <- eta - solve(Hf, S)
      }
      E <- .fa(th, eta, s, times); f <- E$f
      Phi <- 0.5 * sum((y - f)^2 / R0 + log(R0)) + 0.5 * as.numeric(t(eta) %*% Oi %*% eta) + 0.5 * ldOm
      Ht <- Oi; for (l in 1:neta) for (m in 1:neta) Ht[l, m] <- Ht[l, m] + sum((1/R0) * E$a[, l] * E$a[, m])
      tot <- tot + Phi + 0.5 * log(det(Ht))
    }
    tot
  }
  psi0 <- c(as.numeric(fitF$theta[c("tka", "tcl", "tv", "add.sd", "prop.sd")]), diag(Om0))
  np <- length(psi0); h <- pmax(abs(psi0), 0.5) * 5e-5; H <- matrix(0, np, np); f0 <- objFOCE(psi0)
  for (i in 1:np) { ei <- numeric(np); ei[i] <- h[i]
    H[i, i] <- (objFOCE(psi0 + 2*ei) - 2*f0 + objFOCE(psi0 - 2*ei)) / (4 * h[i]^2) }
  for (i in 1:(np-1)) for (j in (i+1):np) { ei <- numeric(np); ei[i] <- h[i]; ej <- numeric(np); ej[j] <- h[j]
    H[i, j] <- H[j, i] <- (objFOCE(psi0+ei+ej) - objFOCE(psi0+ei-ej) - objFOCE(psi0-ei+ej) + objFOCE(psi0-ei-ej)) / (4 * h[i] * h[j]) }
  pn <- c("tka", "tcl", "tv", "add.sd", "prop.sd", "om.eta.ka", "om.eta.cl", "om.eta.v")
  dimnames(H) <- list(pn, pn); Ran <- rF$R[pn, pn]
  # the analytic observed-information R reproduces the gold-FD Hessian: exact on every
  # numerically significant entry (rel < 3e-4 on entries above 1% of the matrix norm;
  # tiny entries carry only central-FD roundoff).
  big <- abs(H) > 0.01 * max(abs(H))
  expect_lt(max(abs(Ran[big] - H[big]) / abs(H[big])), 3e-4)
  # and the whole matrix agrees to central-FD accuracy
  expect_lt(max(abs(Ran - H) / (abs(H) + 1e-6)), 3e-3)
  # finite SEs (the identified directions) match the gold FD
  seA <- suppressWarnings(sqrt(diag(solve(Ran)))); seG <- suppressWarnings(sqrt(diag(solve(H))))
  fin <- is.finite(seA) & is.finite(seG)
  expect_gt(sum(fin), 4L)                              # most directions are identified
  expect_equal(unname(seA[fin]), unname(seG[fin]), tolerance = 5e-3)
})

test_that("foce+ (foce='foce+') combined analytic cov matches the live-R gold FD", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  # foce+ keeps the LIVE conditional variance R = R(f(theta, eta-hat)) in the objective
  # while the inner problem still drops dR/deta (truncated gradient).  Validation is
  # against a fully independent central-FD Hessian of the live-R foce+ objective
  # (Phi(R) + 0.5 log|H~_FOCE|, EBEs re-solved to S = sum(-(y-f)/R . a) + Omega^-1 eta = 0
  # with R live inside the Newton).  Combined error makes foce+ differ from both
  # "nonmem" FOCE (frozen R0) and FOCEI (interaction term).
  theo <- nlmixr2data::theo_sd
  fitP <- suppressMessages(nlmixr(.cov_combined, theo, "focei",
            foceiControl(print = 0L, covMethod = "", interaction = FALSE, foce = "foce+", sigdig = 6)))
  rP <- foceiCovAnalytic(fitP)
  expect_false(is.null(rP)); expect_identical(rP$method, "analytic")

  ui <- fitP$finalUi; neta <- 3L; etav <- paste0("ETA_", 1:neta, "_")
  am <- .foceiAnalyticAugModelDirs(ui, etav)
  thNames <- names(fitP$theta)
  thBase <- setNames(as.numeric(fitP$theta[thNames]), paste0("THETA_", seq_along(thNames), "_"))
  iTh <- match(c("tka", "tcl", "tv", "add.sd", "prop.sd"), thNames)
  Om0 <- fitP$omega
  byId <- split(fitP$dataSav, as.character(fitP$dataSav$ID))
  ids <- fitP$eta$ID; idCode <- as.integer(ids)
  eta0m <- as.matrix(fitP$eta[, c("eta.ka", "eta.cl", "eta.v")])
  subj <- lapply(seq_along(ids), function(i) {
    s <- byId[[as.character(idCode[i])]]; obs <- s[s$EVID == 0, , drop = FALSE]
    list(s = s, times = obs$TIME, y = obs$DV, eta0 = eta0m[i, ]) })
  .fa <- function(th, eta, s, times) .foceiAnalyticSolveFA(am, c(th, setNames(eta, etav)), s, times, tol = 1e-12)
  objFOCEP <- function(psi) {
    th <- thBase; th[iTh] <- psi[1:5]; sa <- psi[4]; sp <- psi[5]
    Om <- Om0; diag(Om) <- psi[6:8]; Oi <- solve(Om); ldOm <- log(det(Om)); tot <- 0
    for (sj in subj) {
      y <- sj$y; s <- sj$s; times <- sj$times
      eta <- sj$eta0
      for (it in 1:100) {
        E <- .fa(th, eta, s, times); if (is.null(E)) return(NA_real_)
        R <- sa^2 + sp^2 * E$f^2                          # live conditional variance
        q0 <- -(y - E$f) / R
        q1 <- 1 / R + (y - E$f) * (2 * sp^2 * E$f) / R^2  # dq0/df with live R
        S <- as.numeric(Oi %*% eta); for (l in 1:neta) S[l] <- S[l] + sum(q0 * E$a[, l])
        if (max(abs(S)) < 1e-12) break
        Hf <- Oi; for (l in 1:neta) for (m in 1:neta) Hf[l, m] <- Hf[l, m] + sum(q1 * E$a[, l] * E$a[, m] + q0 * E$A[, l, m])
        eta <- eta - solve(Hf, S)
      }
      E <- .fa(th, eta, s, times); f <- E$f
      R <- sa^2 + sp^2 * f^2
      Phi <- 0.5 * sum((y - f)^2 / R + log(R)) + 0.5 * as.numeric(t(eta) %*% Oi %*% eta) + 0.5 * ldOm
      Ht <- Oi; for (l in 1:neta) for (m in 1:neta) Ht[l, m] <- Ht[l, m] + sum((1 / R) * E$a[, l] * E$a[, m])
      tot <- tot + Phi + 0.5 * log(det(Ht))
    }
    tot
  }
  psi0 <- c(as.numeric(fitP$theta[c("tka", "tcl", "tv", "add.sd", "prop.sd")]), diag(Om0))
  np <- length(psi0); h <- pmax(abs(psi0), 0.5) * 5e-5; H <- matrix(0, np, np); f0 <- objFOCEP(psi0)
  for (i in 1:np) { ei <- numeric(np); ei[i] <- h[i]
    H[i, i] <- (objFOCEP(psi0 + 2*ei) - 2*f0 + objFOCEP(psi0 - 2*ei)) / (4 * h[i]^2) }
  for (i in 1:(np-1)) for (j in (i+1):np) { ei <- numeric(np); ei[i] <- h[i]; ej <- numeric(np); ej[j] <- h[j]
    H[i, j] <- H[j, i] <- (objFOCEP(psi0+ei+ej) - objFOCEP(psi0+ei-ej) - objFOCEP(psi0-ei+ej) + objFOCEP(psi0-ei-ej)) / (4 * h[i] * h[j]) }
  pn <- c("tka", "tcl", "tv", "add.sd", "prop.sd", "om.eta.ka", "om.eta.cl", "om.eta.v")
  dimnames(H) <- list(pn, pn); Ran <- rP$R[pn, pn]
  # exact on every numerically significant entry (entries below 1% of the matrix norm
  # sit at ~1e-7 of it and carry only central-FD roundoff, so they get a norm-scaled bound)
  big <- abs(H) > 0.01 * max(abs(H))
  expect_lt(max(abs(Ran[big] - H[big]) / abs(H[big])), 3e-4)
  expect_lt(max(abs(Ran - H)), 1e-5 * max(abs(H)))
  # finite SEs (the identified directions) match the gold FD
  seA <- suppressWarnings(sqrt(diag(solve(Ran)))); seG <- suppressWarnings(sqrt(diag(solve(H))))
  fin <- is.finite(seA) & is.finite(seG)
  expect_gt(sum(fin), 4L)                              # most directions are identified
  expect_equal(unname(seA[fin]), unname(seG[fin]), tolerance = 5e-3)
})

test_that("est='focep' installs the full analytic covariance", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")
  fit <- suppressWarnings(suppressMessages(nlmixr(.cov_one_cmt, nlmixr2data::theo_sd, "focep",
            focepControl(print = 0L, covMethod = "analytic", covFull = TRUE))))
  expect_identical(fit$covMethod, "analytic")
  expect_true(any(grepl("^om\\.", rownames(fit$cov))))
  .se <- sqrt(diag(fit$cov))
  expect_true(all(is.finite(.se)) && all(.se > 0))
})
