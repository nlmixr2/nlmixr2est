nmTest({
  # #1132: a residual variance below sqrt(eps) used to be REPLACED by 1, putting a
  # ~+16 cliff per observation in the FOCEi objective where a proportional-error
  # prediction crossed ~1.2e-3.  It is now floored, leaving only a small log|H| step.
  test_that("focei objective is continuous where the variance crosses sqrt(eps)", {
    m <- function() {
      ini({
        lf <- log(0.0012)
        eta.f ~ 0.1
        prop.sd <- 0.1
      })
      model({
        ipred <- exp(lf + eta.f)
        ipred ~ prop(prop.sd)
      })
    }
    d <- data.frame(ID = 1, TIME = c(1, 2), DV = c(0.0012, 0.0011), EVID = 0, AMT = 0)
    sd <- 0.1
    omega <- 0.1
    # eta at which the prediction's variance (sd * f)^2 equals sqrt(eps)
    fCross <- sqrt(sqrt(.Machine$double.eps)) / sd
    etaCross <- log(fCross / 0.0012)

    objAt <- function(eta) {
      ctl <- foceiControl(
        print = 0,
        covMethod = "",
        calcTables = FALSE,
        maxOuterIterations = 0,
        maxInnerIterations = 0,
        etaMat = matrix(eta, 1)
      )
      suppressWarnings(nlmixr2(m, d, "focei", ctl))$objf
    }
    # R-side FOCEi objective: df/deta = f, dR/deta = 2R
    objR <- function(eta) {
      f <- 0.0012 * exp(eta)
      r <- (sd * f)^2
      h <- 1 / omega + nrow(d) * (f^2 / r + 0.5 * (2 * r)^2 / r^2)
      sum(log(r) + (d$DV - f)^2 / r) + eta^2 / omega + log(omega) + log(h)
    }

    # floored side: R is held at sqrt(eps), so dR/deta = 0
    objFloor <- function(eta) {
      f <- 0.0012 * exp(eta)
      r <- sqrt(.Machine$double.eps)
      h <- 1 / omega + nrow(d) * f^2 / r
      sum(log(r) + (d$DV - f)^2 / r) + eta^2 / omega + log(omega) + log(h)
    }

    below <- objAt(etaCross - 1e-4)
    above <- objAt(etaCross + 1e-4)
    # differential pair straddling the threshold: the old code jumped ~32 here;
    # only the small log|H| step from dropping dR/deta remains
    expect_lt(abs(above - below), 0.05)
    # each side matches its closed form
    expect_equal(above, objR(etaCross + 1e-4), tolerance = 1e-4)
    expect_equal(below, objFloor(etaCross - 1e-4), tolerance = 1e-4)
  })

  # The fast=TRUE analytic outer gradient must differentiate the SAME floored
  # objective: it used to differentiate log(R_raw) at a floored observation, and one
  # such subject stopped a Michaelis-Menten fit 63 OFV short of the optimum (#1132).
  test_that("fast=TRUE analytic gradient matches central differences at a floored R", {
    skip_on_cran()
    m <- function() {
      ini({
        lf <- 0
        lk <- log(0.2)
        eta.f ~ 0.1
        eta.k ~ 0.1
        prop.sd <- 0.1
      })
      model({
        ipred <- exp(lf + eta.f) * exp(-exp(lk + eta.k) * TIME)
        ipred ~ prop(prop.sd)
      })
    }
    .testSeed(1132)
    obsT <- c(1, 2, 4, 8, 24, 48)
    d <- do.call(rbind, lapply(1:8, function(i) {
      data.frame(ID = i, TIME = obsT, AMT = 0, EVID = 0,
                 DV = exp(-0.2 * obsT) * exp(rnorm(1, 0, 0.3)) * (1 + 0.1 * rnorm(6)))
    }))
    ctl <- function(fast) {
      foceiControl(print = 0L, covMethod = "", fast = fast, sigdig = 4,
                   maxOuterIterations = 0L, maxInnerIterations = 500L)
    }
    fit <- suppressMessages(suppressWarnings(nlmixr2(m, d, "focei", ctl(TRUE))))
    # the design must actually exercise the floor
    expect_true(any((0.1 * fit$IPRED)^2 < sqrt(.Machine$double.eps)))
    g <- .foceiGradDirect(fit)
    expect_false(is.null(g))
    expect_gt(fit$env$nAnalyticGradDirect, 0)
    base <- fixef(fit)
    ofvAt <- function(nm, val) {
      ui2 <- do.call(rxode2::ini, c(list(fit$finalUi), setNames(list(val), nm)))
      suppressMessages(suppressWarnings(nlmixr2(ui2, d, "focei", ctl(FALSE))))$objf
    }
    fd <- vapply(names(base), function(nm) {
      h <- 1e-4 * max(abs(base[[nm]]), 0.05)
      (ofvAt(nm, base[nm] + h) - ofvAt(nm, base[nm] - h)) / (2 * h)
    }, numeric(1))
    expect_equal(unname(g[names(base)]), unname(fd), tolerance = 0.02)
  })
})
