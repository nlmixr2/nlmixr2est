# General log-likelihood endpoints for est="rpem" (ll(name) ~ <per-obs loglik>), the
# same saemix-style mechanism saem/fsaem use: the model returns the per-observation
# log-likelihood as its prediction, so the E-step uses it directly as the observation
# loss and there is no residual error parameter (errType 7).  Fixed-effect parameters of
# the likelihood expression are estimated by the numeric re-solve (structural-beta) M-step.

.rpemExpTte <- function() {
  ini({ tlam <- log(25); eta.lam ~ 0.2 })
  model({ lam <- exp(tlam + eta.lam); ll(dv) ~ -log(lam) - time / lam })
}

.rpemMkTte <- function(seed = 1L, n = 150L, meanT = 40) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n), function(i) {
    lami <- meanT * exp(rnorm(1, 0, sqrt(0.15)))
    data.frame(ID = i, TIME = lami * -log(runif(1)), DV = 1, EVID = 0, CMT = 1)
  }))
}

test_that("RPEM classifies a general log-likelihood endpoint as errType 7 (no residual)", {
  ui <- rxode2::rxUiDecompress(rxode2::rxode2(.rpemExpTte))
  expect_equal(as.character(ui$predDf$distribution[1]), "LL")
  cl <- .rpemClassify(ui)
  expect_equal(cl$errType, 7L)
  expect_true(is.na(cl$addSdIdx))          # no residual error parameter
  expect_equal(cl$errName, "ll")
})

test_that("RPEM fits a general log-likelihood (exponential TTE) endpoint", {
  skip_on_cran()
  skip_on_ci()  # heavy: multi-iteration RPEM loop

  ui <- rxode2::rxUiDecompress(rxode2::rxode2(.rpemExpTte))
  rf <- .rpemFit(ui, .rpemMkTte(1L),
                 rpemControl(nGauss = 500L, nMH = 60000L, mhBurn = 6000L, niter = 40L,
                             collect = 15L, seed = 1L, cores = 4L))
  # recovers the population mean (true 40) from a poor start (25); no residual sd
  expect_equal(exp(unname(rf$mu["tlam"])), 40, tolerance = 0.2)
  expect_gt(rf$omega[1], 0)
  expect_true(is.na(rf$addSd))
})

test_that("RPEM recovers a fixed-effect likelihood parameter (Weibull shape) as a structural beta", {
  skip_on_cran()
  skip_on_ci()

  weiTte <- function() {
    ini({ tlam <- log(30); lk <- log(1.0); eta.lam ~ 0.2 })
    model({ lam <- exp(tlam + eta.lam); k <- exp(lk)
            ll(dv) ~ lk + (k - 1) * log(time) - k * log(lam) - (time / lam)^k })
  }
  mkWei <- function(seed = 1L, n = 200L, meanT = 40, k = 1.5) {
    set.seed(seed)
    do.call(rbind, lapply(seq_len(n), function(i) {
      lami <- meanT * exp(rnorm(1, 0, sqrt(0.15)))
      data.frame(ID = i, TIME = lami * (-log(runif(1)))^(1 / k), DV = 1, EVID = 0, CMT = 1)
    }))
  }
  ui <- rxode2::rxUiDecompress(rxode2::rxode2(weiTte))
  cl <- .rpemClassify(ui)
  expect_true("lk" %in% cl$thetaNames[cl$structIdx + 1L])   # shape -> structural set
  rf <- .rpemFit(ui, mkWei(1L),
                 rpemControl(nGauss = 500L, nMH = 60000L, mhBurn = 6000L, niter = 40L,
                             collect = 15L, seed = 1L, cores = 4L))
  expect_equal(exp(unname(rf$mu["tlam"])), 40, tolerance = 0.2)
  expect_equal(exp(unname(rf$struct["lk"])), 1.5, tolerance = 0.25)   # Weibull shape
})

test_that("likLbfgs refines a bounded likelihood parameter and clamps to its declared bound", {
  skip_on_cran()
  skip_on_ci()

  # true Weibull shape 4, but the model DECLARES an upper bound of 3: the box-constrained
  # L-BFGS-B refinement (likLbfgs, saem ind.fix10 style) must clamp at 3, whereas the
  # unbounded damped-Newton re-solve overshoots past it.
  weiB <- function() {
    ini({ tlam <- log(38); shape <- c(0.5, 2.0, 3.0); eta.lam ~ 0.2 })
    model({ lam <- exp(tlam + eta.lam)
            ll(dv) ~ log(shape) + (shape - 1) * log(time) - shape * log(lam) - (time / lam)^shape })
  }
  mkWei <- function(seed = 1L, n = 200L, meanT = 40, k = 4.0) {
    set.seed(seed)
    do.call(rbind, lapply(seq_len(n), function(i) {
      lami <- meanT * exp(rnorm(1, 0, sqrt(0.15)))
      data.frame(ID = i, TIME = lami * (-log(runif(1)))^(1 / k), DV = 1, EVID = 0, CMT = 1)
    }))
  }
  ui <- rxode2::rxUiDecompress(rxode2::rxode2(weiB))
  cl <- .rpemClassify(ui)
  expect_equal(cl$structNbd, 2L)            # shape has both bounds
  expect_equal(cl$structUpper, 3)
  d <- mkWei(1L, k = 4.0)
  ctl <- function(lb) rpemControl(nGauss = 400L, nMH = 50000L, mhBurn = 5000L, niter = 35L,
                                  collect = 12L, seed = 1L, cores = 4L, likLbfgs = lb)
  rfC <- .rpemFit(ui, d, ctl(TRUE))
  rfN <- .rpemFit(ui, d, ctl(FALSE))
  # box-constrained refinement never leaves the declared bound
  expect_lte(unname(rfC$struct["shape"]), 3 + 1e-6)
  expect_equal(unname(rfC$struct["shape"]), 3, tolerance = 1e-4)
  # the unbounded re-solve overshoots past the declared bound
  expect_gt(unname(rfN$struct["shape"]), 3)

  # the C++ cLoop applies the same box-constrained refinement, so a bounded ll()
  # parameter also clamps at its declared bound when the whole loop runs in C++
  rfCpp <- .rpemFit(ui, d, rpemControl(nGauss = 400L, nMH = 50000L, mhBurn = 5000L,
                                       niter = 35L, collect = 12L, seed = 1L, cores = 4L,
                                       cLoop = TRUE))
  expect_equal(unname(rfCpp$struct["shape"]), 3, tolerance = 1e-4)
})

test_that("a general log-likelihood endpoint runs entirely in the C++ cLoop", {
  skip_on_cran()
  skip_on_ci()

  ui <- rxode2::rxUiDecompress(rxode2::rxode2(.rpemExpTte))
  d <- .rpemMkTte(1L)
  ctl <- function(cl) rpemControl(nGauss = 400L, nMH = 50000L, mhBurn = 5000L, niter = 30L,
                                  collect = 12L, seed = 1L, cores = 4L, cLoop = cl)
  rfR <- .rpemFit(ui, d, ctl(FALSE))
  rfC <- .rpemFit(ui, d, ctl(TRUE))
  # the C++ loop matches the R loop and recovers the mean; no residual sd
  expect_equal(exp(unname(rfC$mu["tlam"])), 40, tolerance = 0.2)
  expect_equal(unname(rfC$mu["tlam"]), unname(rfR$mu["tlam"]), tolerance = 0.02)
  expect_true(is.na(rfC$addSd))

  # dynamic-iteration stable: a longer C++ run shares the shorter run's per-iteration prefix
  cl <- .rpemClassify(ui)
  .nm <- c(paste0("THETA[", seq_len(cl$nTheta), "]"), paste0("ETA[", seq_len(cl$nEta), "]"))
  e <- new.env(); e$predOnly <- ui$rpemRxModel$predOnly
  e$rxControl <- rxode2::rxControl(atol = 1e-8, rtol = 1e-8, cores = 2L)
  e$param <- stats::setNames(cl$base, .nm); e$data <- d
  runN <- function(ni) rpemEMLoopK1(e, cl$base, cl$etaIdx, cl$muIdx, -1L, 7L, cl$mu0,
    diag(as.matrix(cl$omega0)), 0.0, c(-1L, -1L, -1L), c(0, 0, 0), integer(0), numeric(0),
    ni, 300L, 2L, 40000L, 4000L, 7L, matrix(0, 0, 0), integer(0),
    numeric(0), numeric(0), integer(0), 0L, 12L, 5L, 1e7, 0, 20L, 1.0)
  a <- runN(15L); b <- runN(25L)
  expect_equal(a$muTrace[1:15, , drop = FALSE], b$muTrace[1:15, , drop = FALSE])
  expect_equal(as.numeric(a$lnL)[1:15], as.numeric(b$lnL)[1:15])
})
