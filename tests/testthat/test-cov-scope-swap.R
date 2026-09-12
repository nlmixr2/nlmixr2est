# Both covariance shapes (foceiControl(covFull=)) are named and cached, so
# setCov() swaps between the theta-only and the full theta+sigma+Omega matrix
# without recomputing either.  The invariant every swap must keep is
# parFixedDf$SE == sqrt(diag(fit$cov)).

nmTest({
  .seMatchesCov <- function(fit) {
    .se <- fit$parFixedDf$SE
    names(.se) <- rownames(fit$parFixedDf)
    .d <- sqrt(diag(fit$cov))
    .n <- intersect(names(.se), names(.d))
    expect_true(length(.n) > 0L)
    expect_equal(unname(.se[.n]), unname(.d[.n]))
  }

  .oneCmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1.0
      tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("a covFull focei fit names the full shape and caches the theta-only one", {
    .fit <- suppressWarnings(nlmixr2(
      .oneCmt, nlmixr2data::theo_sd, est = "focei",
      control = foceiControl(print = 0, covMethod = "s", calcTables = FALSE)))
    expect_identical(.fit$covMethod, "s (full)")
    expect_true(nrow(.fit$cov) > nrow(.fit$parFixedDf))   # theta + Omega
    expect_true("s" %in% names(.fit$env$covList))
    .seMatchesCov(.fit)

    .cached <- .fit$env$covList[["s"]]
    .fullSe <- sqrt(diag(.fit$cov))

    # swapping to the theta-only shape reinstalls the CACHED matrix (not a recompute)
    setCov(.fit, "s")
    expect_identical(.fit$covMethod, "s")
    expect_equal(unname(.fit$cov), unname(.cached))
    .seMatchesCov(.fit)
    # the FD shapes are different estimators: solve(S_theta) vs the theta block of
    # solve(S_full), which also carries the Omega estimation uncertainty
    .thetaSe <- sqrt(diag(.fit$cov))
    expect_false(isTRUE(all.equal(unname(.thetaSe),
                                  unname(.fullSe[names(.thetaSe)]))))

    # ... and the full shape is now the cached one, so the round trip is exact
    expect_true("s (full)" %in% names(.fit$env$covList))
    setCov(.fit, "s (full)")
    expect_identical(.fit$covMethod, "s (full)")
    expect_equal(unname(sqrt(diag(.fit$cov))), unname(.fullSe))
    .seMatchesCov(.fit)
  })

  test_that("setCov() refuses to switch to the shape already installed", {
    .fit <- suppressWarnings(nlmixr2(
      .oneCmt, nlmixr2data::theo_sd, est = "focei",
      control = foceiControl(print = 0, covMethod = "s", calcTables = FALSE)))
    expect_error(setCov(.fit, "s (full)"), "no need to switch")
    expect_error(setCov(.fit, "nonesuch"), "have not been calculated")
  })

  test_that("covFull=FALSE reports the unqualified name and caches nothing", {
    .fit <- suppressWarnings(nlmixr2(
      .oneCmt, nlmixr2data::theo_sd, est = "focei",
      control = foceiControl(print = 0, covMethod = "s", covFull = FALSE,
                             calcTables = FALSE)))
    expect_identical(.fit$covMethod, "s")
    expect_equal(nrow(.fit$cov), nrow(.fit$parFixedDf))
    # only one shape was computed, so there is nothing to swap to
    expect_false(exists("covList", envir = .fit$env, inherits = FALSE))
    .seMatchesCov(.fit)
  })

  .odeCmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1.0
      tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  test_that("the analytic shapes agree on the theta SEs and swap both ways", {
    .fit <- suppressWarnings(nlmixr2(
      .odeCmt, nlmixr2data::theo_sd, est = "focei",
      control = foceiControl(print = 0, calcTables = FALSE)))

    setCov(.fit, "analytic")
    expect_identical(.fit$covMethod, "analytic")
    expect_equal(nrow(.fit$cov), nrow(.fit$parFixedDf))
    .thetaSe <- sqrt(diag(.fit$cov))
    .seMatchesCov(.fit)
    # assembling the analytic covariance produced both shapes; the other is cached
    expect_true("analytic (full)" %in% names(.fit$env$covList))
    .cachedFull <- .fit$env$covList[["analytic (full)"]]

    setCov(.fit, "analytic (full)")
    expect_identical(.fit$covMethod, "analytic (full)")
    expect_equal(unname(.fit$cov), unname(.cachedFull))
    .seMatchesCov(.fit)
    # the analytic assembly is always full, so the theta-only shape is a submatrix
    # of the inverse -- unlike the FD path, the theta SEs agree
    .fullSe <- sqrt(diag(.fit$cov))
    expect_equal(unname(.thetaSe), unname(.fullSe[names(.thetaSe)]))
  })
})
