# The finite-difference covariance options foceiControl() and rsControl() share.

test_that("foceiControl() and rsControl() take the same finite-difference covariance options", {
  # gillStepCov is the factor the Gill search grows its step by
  expect_error(foceiControl(gillStepCov = 0.5), "gillStepCov")
  expect_error(rsControl(gillStepCov = 0.5), "gillStepCov")
  expect_identical(foceiControl(gillStepCov = 1)$gillStepCov, 1)
  expect_identical(rsControl(gillStepCov = 1)$gillStepCov, 1)
  # covSmall is a single finite threshold
  expect_error(foceiControl(covSmall = Inf), "covSmall")
  expect_error(rsControl(covSmall = Inf), "covSmall")
  expect_error(foceiControl(covSmall = c(1e-5, 1e-6)), "covSmall")
  expect_error(rsControl(covSmall = c(1e-5, 1e-6)), "covSmall")
  # the flags are logical or 0/1, which foceiControl() stores
  expect_identical(
    unclass(rsControl(covGillF = 0L, rmatNorm = 1, smatNorm = TRUE)),
    list(covGillF = 0L, rmatNorm = 1, smatNorm = TRUE)
  )
  expect_identical(
    foceiControl(covGillF = 0L, rmatNorm = 1, smatNorm = TRUE)[c("covGillF", "rmatNorm", "smatNorm")],
    list(covGillF = 0L, rmatNorm = 1L, smatNorm = 1L)
  )
  expect_error(rsControl(rmatNorm = 2), "rmatNorm")
  expect_error(foceiControl(rmatNorm = 2), "rmatNorm")
})

test_that("every Gill step factor is a finite number of at least 1", {
  # gillStep (the fit's gradient search), gillStepCov and gillStepCovLlik (the
  # covariance step) all grow the Gill (1983) step by multiplying by the factor
  # and shrink it by dividing
  for (.n in c("gillStep", "gillStepCov", "gillStepCovLlik")) {
    .below <- paste0("Assertion on '", .n, "' failed: Element 1 is not >= 1.")
    .inf <- paste0("Assertion on '", .n, "' failed: Must be finite.")
    expect_error(do.call(foceiControl, setNames(list(0.5), .n)), .below, fixed = TRUE)
    expect_error(do.call(foceiControl, setNames(list(0), .n)), .below, fixed = TRUE)
    expect_error(do.call(foceiControl, setNames(list(Inf), .n)), .inf, fixed = TRUE)
    expect_identical(do.call(foceiControl, setNames(list(1), .n))[[.n]], 1)
  }
  expect_error(rsControl(gillStepCov = Inf), "Assertion on 'gillStepCov' failed: Must be finite.", fixed = TRUE)
  # nlmixr2Gill83() makes the checks foceiControl() makes
  .f <- function(x) sum(x^2)
  expect_error(
    nlmixr2Gill83(.f, c(1, 2), gillStep = 0.5),
    "Assertion on 'gillStep' failed: Element 1 is not >= 1.",
    fixed = TRUE
  )
  expect_error(
    nlmixr2Gill83(.f, c(1, 2), gillStep = Inf),
    "Assertion on 'gillStep' failed: Must be finite.",
    fixed = TRUE
  )
})

test_that("rsControl(covFallback=) is the request's own list, never the fit's", {
  expect_null(rsControl()$covFallback)
  expect_identical(rsControl(covFallback = list(r = "s"))$covFallback, list(r = "s"))
  expect_error(rsControl(covFallback = list(r = "vi")), "cannot fall back to \"vi\"")
  .fit <- new.env(parent = emptyenv())
  .fit$foceiControl <- foceiControl()
  # the fallbacks decide what is installed, not how a covariance is computed, so a
  # cached covariance is reused whatever the list
  expect_identical(setCovOptions(rsControl(covFallback = list(r = "s")), .fit), setCovOptions(rsControl(), .fit))
  .d <- rxode2::rxUiDeparse(rsControl(covFallback = list(r = "s")), "ctl")
  expect_identical(eval(.d[[3]])$covFallback, list(r = "s"))
})

nmTest({
  test_that("setCov() falls back only as rsControl(covFallback=) lists", {
    skip_on_cran()
    .m <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    # an "s" fit (no R computed, so no "r" to swap to) and an "r" request whose refit
    # can only give "s"
    .fit <- .nlmixr(.m, theo_sd, "focei", foceiControl(print = 0, calcTables = FALSE, covMethod = "s", covFull = FALSE))
    expect_null(.fit$env$covList$r)
    .s <- suppressMessages(.setCovRefit(.fit, covMethod = "s", covFull = FALSE))
    expect_identical(.s$covMethod, "s")
    .args <- new.env(parent = emptyenv())
    local_mocked_bindings(.setCovRefit = function(obj, ...) {
      .args$fallback <- list(...)$covFallback
      .s
    })
    expect_error(suppressMessages(setCov(.fit, "r")), "\"r\" could not be computed")
    # the refit was given no fallback
    expect_identical(.args$fallback, list())
    expect_warning(
      suppressMessages(setCov(.fit, "r", control = rsControl(covFallback = list(r = "s")))),
      "\"s\" covariance installed instead of the requested \"r\""
    )
    expect_identical(.args$fallback, list(r = "s"))
    expect_identical(.fit$covMethod, "s")
    expect_identical(.fit$cov, .s$cov)
  })
})

test_that("a SAEM fit's \"sa\" recompute continues its chains with no warm-up iterations", {
  .st <- list(phiM = matrix(seq_len(18) / 10, 6, 3), sigma2 = 0.5)
  .warm <- .covEngineControl("sa", saControl(nBurn = 7L, nEm = 8L), .st)
  expect_identical(.warm$mcmc$niter, c(0L, 0L))
  expect_identical(.warm$saemWarmState, .st)
  expect_true(.warm$saemHoldPar)
  .cold <- .covEngineControl("sa", saControl(nBurn = 7L, nEm = 8L))
  expect_identical(.cold$mcmc$niter, c(7L, 8L))
  expect_null(.cold$saemWarmState)
  expect_error(saControl(warmStart = NA), "warmStart")
  expect_error(saemControl(saemWarmState = list(phiM = matrix(NA_real_, 2, 2))), "phiM")
  expect_error(saemControl(saemWarmState = list(phiM = .st$phiM, sigma2 = -1)), "sigma2")
  expect_error(saemControl(saemWarmState = list(phiM = .st$phiM, mpostPhi = matrix(NA_real_, 2, 3))), "mpostPhi")
})

test_that(".saemWarmCfg() installs a chain state of the right shape, its statistics and the residual statistic", {
  # 2 subjects x 3 chains, phi columns 1:2 mu-referenced (i1) and 3 not (i0);
  # endpoints: additive (4 obs), proportional (3 obs), combined (5 obs)
  .cfg <- list(
    phiM = matrix(0, 6, 3),
    i1 = 0:1,
    i0 = 2L,
    N = 2L,
    nmc = 3L,
    nMix = 1L,
    res.mod = c(1, 2, 4),
    ares = c(10, 0, 10),
    bres = c(0, 1, 1),
    y_offset = c(0, 4, 7, 12),
    res_offset = c(0L, 1L, 2L, 4L),
    resValue = c(0.5, 0.1, 0.2, 0.3),
    Gamma2_phi1 = diag(c(0.4, 0.3)),
    Gamma2_phi1fixedIx = matrix(1L, 2, 2),
    Gamma2_phi1fixedValues = matrix(c(0.5, 0.1, 0.1, 0.2), 2, 2)
  )
  .ph <- matrix(c(1, 2, 3, 4, 5, 6, 10, 20, 30, 40, 50, 60, 7, 8, 9, 10, 11, 12), 6, 3)
  .w <- .saemWarmCfg(.cfg, list(phiM = .ph))
  expect_identical(.w$phiM, .ph)
  # row i + k * N is subject i of chain k: subject 1 holds rows 1, 3, 5
  expect_equal(.w$statphi11, matrix(c(3, 4, 30, 40), 2, 2))
  expect_equal(.w$statphi01, matrix(c(9, 10), 2, 1))
  # the kernel's second moments are averaged over the chains (Statphi12 / nmc)
  expect_equal(.w$statphi12, crossprod(.ph[, 1:2]) / 3)
  expect_equal(.w$statphi02, crossprod(.ph[, 3, drop = FALSE]) / 3)
  # the residual parameters start at their held values, not the placeholder 10
  expect_equal(.w$ares, c(0.5, 0, 0.2))
  expect_equal(.w$bres, c(0, 0.1, 0.3))
  # without the fit's sigma2, statrese / n is the held variance
  expect_equal(.w$statrese, c(4 * 0.25, 3 * 0.01, 5))
  # with it, statrese / n is the sigma2 the fit's Louis residual score last read
  .s <- .saemWarmCfg(.cfg, list(phiM = .ph, sigma2 = c(0.3, 0.02, 0.9)))
  expect_equal(.s$statrese, c(4 * 0.3, 3 * 0.02, 5 * 0.9))
  expect_equal(.saemWarmCfg(.cfg, list(phiM = .ph, sigma2 = 1))$statrese, .w$statrese)
  # Omega starts at the held values, covariances included
  expect_equal(.w$Gamma2_phi1, matrix(c(0.5, 0.1, 0.1, 0.2), 2, 2))
  # the posterior means the linearized FIM is taken at; none without them
  expect_null(.w$mpost_phi)
  .mp <- matrix(c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), 2, 3)
  expect_identical(.saemWarmCfg(.cfg, list(phiM = .ph, mpostPhi = .mp))$mpost_phi, .mp)
  expect_null(.saemWarmCfg(.cfg, list(phiM = .ph, mpostPhi = .mp[, 1:2]))$mpost_phi)
  expect_identical(.saemWarmCfg(.cfg, NULL), .cfg)
  expect_identical(.saemWarmCfg(.cfg, list(phiM = .ph[1:4, ])), .cfg)
  .mix <- modifyList(.cfg, list(nMix = 2L))
  expect_identical(.saemWarmCfg(.mix, list(phiM = .ph)), .mix)
})

test_that(".saemHoldCfg() gives the covariance phase an estimated sigma's held value", {
  # endpoints: additive estimated (4 obs), proportional fixed (3 obs), combined (5 obs),
  # additive estimated (2 obs)
  .cfg <- list(
    nlambda1 = 2L,
    nlambda0 = 1L,
    covstruct1 = diag(2),
    res.mod = c(1, 2, 4, 1),
    y_offset = c(0, 4, 7, 12, 14),
    res_offset = c(0L, 1L, 2L, 4L),
    resFixed = c(0L, 1L, 0L, 0L, 0L),
    resValue = c(0.5, 0.1, 0.2, 0.3, 0.7)
  )
  .h <- .saemHoldCfg(.cfg)
  expect_equal(.h$statreseCov, c(4 * 0.25, NA, NA, 2 * 0.49))
  expect_identical(.h$resFixed, rep(1L, 5))
  expect_identical(.h$resKeep, c(0L, 2L, 3L, 4L))
})
