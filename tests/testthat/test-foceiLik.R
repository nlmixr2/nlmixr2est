# General FOCE-family per-subject log-likelihood (issue #414): the
# foceiLikLoad/foceiLikRun/foceiLikUnload lifecycle.  Always-run core test: one
# small analytic (non-ODE) model so it compiles fast and the prediction is
# hand-computable.

.foceiLikMod <- function() {
  ini({
    tcl <- 1
    tv <- 3
    add.sd <- 0.5
    eta.cl ~ 0.1
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    cp <- 100 / v * exp(-cl / v * time)
    cp ~ add(add.sd)
  })
}

.foceiLikData <- function() {
  set.seed(42)
  do.call(rbind, lapply(1:4, function(id) {
    tt <- c(0.5, 1, 2, 4, 8)
    data.frame(ID = id, TIME = tt,
               DV = 5 * exp(-0.05 * tt) + stats::rnorm(length(tt), 0, 0.5),
               AMT = 0, EVID = 0)
  }))
}

test_that("foceiLikLoad() loads, reports dimensions, and guards double loading", {
  skip_on_cran()
  d <- .foceiLikData()
  h <- foceiLikLoad(.foceiLikMod, d, "focei")
  on.exit(foceiLikUnload(), add = TRUE)
  expect_equal(h$nid, 4L)
  expect_equal(h$neta, 1L)
  expect_equal(h$npars, 4L)   # tcl, tv, add.sd + one omega element
  expect_equal(h$etaNames, "eta.cl")
  expect_equal(h$thetaNames, c("tcl", "tv", "add.sd"))
  expect_equal(h$likelihood, "focei")
  # initPar is the estimation-scale vector at the ini() estimates; the omega
  # element uses the "sqrt" diag.xform (omega^(-1/4) for a diagonal omega)
  expect_equal(h$initPar[1:3], c(1, 3, 0.5))
  expect_equal(h$initPar[4], 0.1^(-0.25))
  # loading a second system while one is loaded must error
  expect_error(foceiLikLoad(.foceiLikMod, d, "focei"), "already loaded")
})

test_that("foceiLikUnload() frees and allows a reload; run errors when unloaded", {
  skip_on_cran()
  d <- .foceiLikData()
  h <- foceiLikLoad(.foceiLikMod, d, "focei")
  expect_true(foceiLikUnload())
  # nothing loaded -> unload is a no-op, run errors
  expect_false(foceiLikUnload())
  expect_error(foceiLikRun(h$initPar, matrix(0, 4, 1)), "no general likelihood system loaded")
  # reload works after an unload
  h2 <- foceiLikLoad(.foceiLikMod, d, "foce")
  on.exit(foceiLikUnload(), add = TRUE)
  expect_equal(h2$npars, 4L)
})

test_that("foceiLikRun() returns per-id log-likelihoods on the right scale", {
  skip_on_cran()
  d <- .foceiLikData()
  h <- foceiLikLoad(.foceiLikMod, d, "focei")
  on.exit(foceiLikUnload(), add = TRUE)
  eta0 <- matrix(0, h$nid, h$neta)
  llData <- foceiLikRun(h$initPar, eta0, type = "cond")
  expect_length(llData, 4L)
  expect_true(all(is.finite(llData)))
  expect_equal(names(llData), as.character(h$idLvl))
  # "joint" is the default type, and it differs from "cond" (by the eta prior)
  expect_equal(foceiLikRun(h$initPar, eta0),
               foceiLikRun(h$initPar, eta0, type = "joint"))
  expect_false(isTRUE(all.equal(as.numeric(llData),
                               as.numeric(foceiLikRun(h$initPar, eta0)))))
  # Reference at eta=0 in nlmixr2's residual convention:
  # -0.5*err^2/r - 0.5*log(r) per observation (Gaussian 2*pi constant omitted).
  # The remaining difference is the inner (sensitivity) model solve vs this
  # closed-form evaluation, hence the loose tolerance.
  th <- h$initPar
  ref <- vapply(1:4, function(id) {
    di <- d[d$ID == id, ]
    cp <- 100 / exp(th[2]) * exp(-exp(th[1]) / exp(th[2]) * di$TIME)
    sum(-0.5 * (di$DV - cp)^2 / th[3]^2 - 0.5 * log(th[3]^2))
  }, numeric(1))
  expect_equal(as.numeric(llData), ref, tolerance = 1e-2)
})

test_that("foceiLikRun(type='joint') adds exactly the Gaussian eta prior", {
  skip_on_cran()
  d <- .foceiLikData()
  h <- foceiLikLoad(.foceiLikMod, d, "focei")
  on.exit(foceiLikUnload(), add = TRUE)
  omega <- 0.1
  for (etaVal in list(rep(0, 4), c(-0.3, 0.1, 0.25, 0.4))) {
    eta <- matrix(etaVal, h$nid, h$neta)
    llD <- foceiLikRun(h$initPar, eta, type = "cond")
    llO <- foceiLikRun(h$initPar, eta, type = "joint")
    # joint - cond must be the fully normalized N(0, Omega) log density.  The
    # prior is built from the engine's own omegaInv/logDetOmegaInv5 (the same
    # Omega the inner likelihood uses), whose representation differs from the
    # nominal ini() value by ~4e-5 relative -- hence the tolerance here.
    prior <- stats::dnorm(as.numeric(eta), 0, sqrt(omega), log = TRUE)
    expect_equal(as.numeric(llO - llD), prior, tolerance = 1e-4)
  }
})

test_that("foceiLikRun() is invariant to the thread count", {
  skip_on_cran()
  d <- .foceiLikData()
  h <- foceiLikLoad(.foceiLikMod, d, "focei")
  on.exit(foceiLikUnload(), add = TRUE)
  eta <- matrix(c(-0.2, 0.05, 0.3, -0.1), h$nid, h$neta)
  llS <- foceiLikRun(h$initPar, eta, type = "joint", cores = 1L)
  llP <- foceiLikRun(h$initPar, eta, type = "joint", cores = 4L)
  expect_equal(llS, llP, tolerance = 1e-12)
})

test_that("each likelihood type runs; interaction only matters away from eta=0", {
  skip_on_cran()
  d <- .foceiLikData()
  eta0 <- matrix(0, 4, 1)
  etaNz <- matrix(c(-0.4, 0.2, 0.35, -0.15), 4, 1)
  res <- lapply(c("focei", "focep", "foce"), function(lik) {
    h <- foceiLikLoad(.foceiLikMod, d, lik)
    on.exit(foceiLikUnload(), add = TRUE)
    expect_equal(h$likelihood, lik)
    list(zero = foceiLikRun(h$initPar, eta0, type = "joint"),
         nz = foceiLikRun(h$initPar, etaNz, type = "joint"))
  })
  names(res) <- c("focei", "focep", "foce")
  for (r in res) {
    expect_true(all(is.finite(r$zero)))
    expect_true(all(is.finite(r$nz)))
  }
  # This model has an additive (eta-independent) residual variance, so the
  # interaction/foce variants agree; the likelihood value itself does not
  # depend on how R is evaluated here.
  expect_equal(res$focei$nz, res$foce$nz, tolerance = 1e-8)
  expect_equal(res$focei$nz, res$focep$nz, tolerance = 1e-8)
})

test_that("foceiLikRun() validates theta/eta dimensions", {
  skip_on_cran()
  d <- .foceiLikData()
  h <- foceiLikLoad(.foceiLikMod, d, "focei")
  on.exit(foceiLikUnload(), add = TRUE)
  expect_error(foceiLikRun(h$initPar[-1], matrix(0, 4, 1)), "must have length 4")
  expect_error(foceiLikRun(h$initPar, matrix(0, 4, 2)), "must have 1 columns")
  expect_error(foceiLikRun(h$initPar, matrix(0, 3, 1)), "must have 4 rows")
})

test_that("foceiLikRun() responds to theta changes", {
  skip_on_cran()
  d <- .foceiLikData()
  h <- foceiLikLoad(.foceiLikMod, d, "focei")
  on.exit(foceiLikUnload(), add = TRUE)
  eta0 <- matrix(0, h$nid, h$neta)
  ll0 <- foceiLikRun(h$initPar, eta0, type = "cond")
  # doubling the additive SD must change the data log-likelihood
  th2 <- h$initPar; th2[3] <- th2[3] * 2
  ll2 <- foceiLikRun(th2, eta0, type = "cond")
  expect_false(isTRUE(all.equal(as.numeric(ll0), as.numeric(ll2))))
  # and the theta must actually be re-applied (not cached): returning to the
  # original theta reproduces the original value
  ll0b <- foceiLikRun(h$initPar, eta0, type = "cond")
  expect_equal(ll0, ll0b, tolerance = 1e-12)
})
