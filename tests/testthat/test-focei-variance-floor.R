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

  # A proportional-error design whose late predictions floor the variance.
  .floorMod <- function() {
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
  .floorData <- function() {
    .testSeed(1132)
    obsT <- c(1, 2, 4, 8, 24, 48)
    do.call(
      rbind,
      lapply(1:8, function(i) {
        data.frame(
          ID = i,
          TIME = obsT,
          AMT = 0,
          EVID = 0,
          DV = exp(-0.2 * obsT) * exp(rnorm(1, 0, 0.3)) * (1 + 0.1 * rnorm(6))
        )
      })
    )
  }
  .floorCtl <- function(fast = FALSE, sigdig = 4, maxOuterIterations = 0L) {
    foceiControl(
      print = 0L,
      covMethod = "",
      fast = fast,
      sigdig = sigdig,
      maxOuterIterations = maxOuterIterations,
      maxInnerIterations = 500L
    )
  }
  # objective at fit's estimates with some thetas moved, ETAs re-optimized
  .floorOfv <- function(fit, d, v, sigdig = 4) {
    ui2 <- do.call(rxode2::ini, c(list(fit$finalUi), as.list(v)))
    suppressMessages(suppressWarnings(nlmixr2(ui2, d, "focei", .floorCtl(sigdig = sigdig))))$objf
  }

  # The fast=TRUE analytic outer gradient must differentiate the SAME floored
  # objective: it used to differentiate log(R_raw) at a floored observation, and one
  # such subject stopped a Michaelis-Menten fit 63 OFV short of the optimum (#1132).
  test_that("fast=TRUE analytic gradient matches central differences at a floored R", {
    skip_on_cran()
    d <- .floorData()
    fit <- suppressMessages(suppressWarnings(nlmixr2(.floorMod, d, "focei", .floorCtl(TRUE))))
    # the design must actually exercise the floor
    expect_true(any((0.1 * fit$IPRED)^2 < sqrt(.Machine$double.eps)))
    g <- .foceiGradDirect(fit)
    expect_false(is.null(g))
    expect_gt(fit$env$nAnalyticGradDirect, 0)
    base <- fixef(fit)
    fd <- vapply(
      names(base),
      function(nm) {
        h <- 1e-4 * max(abs(base[[nm]]), 0.05)
        (.floorOfv(fit, d, base[nm] + h) - .floorOfv(fit, d, base[nm] - h)) / (2 * h)
      },
      numeric(1)
    )
    expect_equal(unname(g[names(base)]), unname(fd), tolerance = 0.02)
  })

  # The add/prop covariance assembler writes R as a function of f and cannot express
  # the floor; a floored fit is rerouted to the (f,R) assembler, which reads it.  Its
  # observed information must match the Hessian of the floored objective (it was 28%
  # off in the theta block before).
  test_that("analytic covariance matches the objective's Hessian at a floored R", {
    skip_on_cran()
    d <- .floorData()
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(.floorMod, d, "focei", .floorCtl(sigdig = 6, maxOuterIterations = 1000L))
    ))
    base <- fixef(fit)
    expect_true(any((base[["prop.sd"]] * fit$IPRED)^2 < sqrt(.Machine$double.eps)))
    nm <- names(base)
    .rfr <- new.env()
    .rfr$n <- 0L
    trace(
      ".foceiAnalyticAssembleRFR",
      bquote(assign("n", get("n", envir = .(.rfr)) + 1L, envir = .(.rfr))),
      print = FALSE,
      where = environment(.foceiAnalyticAssembleRFR)
    )
    on.exit(untrace(".foceiAnalyticAssembleRFR", where = environment(.foceiAnalyticAssembleRFR)), add = TRUE)
    r <- suppressWarnings(foceiCovAnalytic(fit))
    expect_gt(.rfr$n, 0L)
    expect_identical(r$method, "analytic")
    h <- 2e-3 * pmax(abs(base), 0.05)
    f0 <- .floorOfv(fit, d, base, sigdig = 6)
    H <- matrix(0, 3, 3)
    for (i in 1:3) {
      for (j in i:3) {
        at <- function(si, sj) {
          v <- base
          v[i] <- v[i] + si * h[i]
          v[j] <- v[j] + sj * h[j]
          .floorOfv(fit, d, v, sigdig = 6)
        }
        H[i, j] <- H[j, i] <- if (i == j) {
          (at(1, 0) - 2 * f0 + at(-1, 0)) / h[i]^2
        } else {
          (at(1, 1) - at(1, -1) - at(-1, 1) + at(-1, -1)) / (4 * h[i] * h[j])
        }
      }
    }
    expect_equal(unname(r$R[nm, nm]), H / 2, tolerance = 1e-3)
  })

  # AGQ quadrature nodes can floor R where eta-hat does not; the (f,R) assembler has
  # no nodes, so such a fit falls back to finite differences.
  test_that("AGQ analytic covariance falls back when only a node floors R", {
    skip_on_cran()
    .testSeed(1132)
    obsT <- c(1, 2, 4, 8, 16, 31)
    d <- do.call(
      rbind,
      lapply(1:8, function(i) {
        data.frame(
          ID = i,
          TIME = obsT,
          AMT = 0,
          EVID = 0,
          DV = exp(-0.2 * obsT) * exp(rnorm(1, 0, 0.3)) * (1 + 0.1 * rnorm(6))
        )
      })
    )
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(.floorMod, d, "agq", foceiControl(print = 0L, covMethod = "", nAGQ = 3, maxOuterIterations = 0L))
    ))
    # eta-hat stays above the floor, so only the node check can see it
    expect_true(all((0.1 * fit$IPRED)^2 > sqrt(.Machine$double.eps)))
    expect_message(r <- foceiCovAnalytic(fit), "floored residual variance")
    expect_null(r)
  })

  # An exact-zero prediction keeps the legacy R = 1; its second derivatives
  # (2 * sp^2 * a * a') are not zero, so the gradient must drop them as well.
  test_that("fast=TRUE analytic gradient matches central differences at R = 0", {
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
        ipred <- exp(lf + eta.f) * TIME * exp(-exp(lk + eta.k) * TIME)
        ipred ~ prop(prop.sd)
      })
    }
    .testSeed(1132)
    obsT <- c(0, 1, 2, 4, 8, 12)
    d <- do.call(
      rbind,
      lapply(1:8, function(i) {
        data.frame(
          ID = i,
          TIME = obsT,
          AMT = 0,
          EVID = 0,
          DV = obsT * exp(-0.2 * obsT) * exp(rnorm(1, 0, 0.3)) * (1 + 0.1 * rnorm(6))
        )
      })
    )
    fit <- suppressMessages(suppressWarnings(nlmixr2(m, d, "focei", .floorCtl(TRUE))))
    expect_true(any(fit$IPRED == 0))
    g <- .foceiGradDirect(fit)
    expect_false(is.null(g))
    expect_gt(fit$env$nAnalyticGradDirect, 0)
    base <- fixef(fit)
    fd <- vapply(
      names(base),
      function(nm) {
        h <- 1e-4 * max(abs(base[[nm]]), 0.05)
        (.floorOfv(fit, d, base[nm] + h) - .floorOfv(fit, d, base[nm] - h)) / (2 * h)
      },
      numeric(1)
    )
    expect_equal(unname(g[names(base)]), unname(fd), tolerance = 0.02)
  })
})
