# Full C++ E-M loop (rpemControl(cLoop = TRUE), design/rpem/12 M5): the additive /
# proportional, diagonal-omega, mu-referenced core runs the whole E-M loop in C++,
# drawing etas from rxode2's per-thread threefry engine with a deterministic
# per-(iteration, subject) seed.  It must recover the same estimates as the R loop,
# be reproducible run-to-run, and produce a full nlmixr2FitData.

test_that("est=rpem cLoop recovers, matches the R loop, and is reproducible", {
  skip_on_cran()
  skip_on_ci()  # heavy: two RPEM fits

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(8); nsub <- 40L; etas <- rnorm(nsub, 0, sqrt(0.3)); obsT <- seq(0.5, 24, by = 1.5)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = etas[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, 0.1)
    d
  }))
  mod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.2; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  ui <- rxode2::rxode2(mod)
  ctlR <- rpemControl(nGauss = 300L, nMH = 60000L, mhBurn = 6000L, niter = 25L,
                      collect = 10L, seed = 42L, cores = 4L, cLoop = FALSE)
  ctlC <- rpemControl(nGauss = 300L, nMH = 60000L, mhBurn = 6000L, niter = 25L,
                      collect = 10L, seed = 42L, cores = 4L, cLoop = TRUE)

  fitC <- suppressMessages(nlmixr2(mod, dat, est = "rpem", control = ctlC))
  expect_s3_class(fitC, "nlmixr2FitData")
  # recovers the truth (tka 0.45, add.sd 0.1)
  expect_equal(unname(fitC$parFixedDf["tka", "Estimate"]), 0.45, tolerance = 0.15)
  expect_equal(unname(fitC$parFixedDf["add.sd", "Estimate"]), 0.1, tolerance = 0.1)

  # matches the R loop closely (both target the same posterior; small Monte-Carlo gap)
  rfR <- .rpemFit(ui, dat, ctlR)
  rfC <- .rpemFit(ui, dat, ctlC)
  expect_equal(unname(rfC$mu["tka"]), unname(rfR$mu["tka"]), tolerance = 0.06)
  expect_equal(rfC$addSd, rfR$addSd, tolerance = 0.06)

  # reproducible run-to-run for a fixed core count (deterministic threefry seeds)
  rfC2 <- .rpemFit(ui, dat, ctlC)
  expect_equal(rfC$mu, rfC2$mu)
  expect_equal(rfC$omega, rfC2$omega)
  expect_equal(rfC$addSd, rfC2$addSd)
})

test_that("est=rpem cLoop matches the R loop for combined / power / TBS residuals", {
  skip_on_cran()
  skip_on_ci()  # heavy: several RPEM fits

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  mkdat <- function(seed, errfun) {
    set.seed(seed); nsub <- 50L; et0 <- rnorm(nsub, 0, sqrt(0.3)); obsT <- seq(0.5, 24, by = 1.5)
    do.call(rbind, lapply(seq_len(nsub), function(i) {
      ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
      s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = et0[i]),
                           events = ev, returnType = "data.frame", addDosing = FALSE)
      d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0; d$DV[o] <- errfun(s$cp)
      d
    }))
  }
  ctlR <- rpemControl(nGauss = 300L, nMH = 60000L, mhBurn = 6000L, niter = 25L,
                      collect = 10L, seed = 42L, cores = 4L, cLoop = FALSE)
  ctlC <- rpemControl(nGauss = 300L, nMH = 60000L, mhBurn = 6000L, niter = 25L,
                      collect = 10L, seed = 42L, cores = 4L, cLoop = TRUE)

  # combined (add + prop): the second residual parameter is re-optimized in C++
  datC <- mkdat(21, function(cp) cp + rnorm(length(cp), 0, sqrt(0.05^2 + (0.1 * cp)^2)))
  modC <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.1; prop.sd <- 0.15; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt()
            cp ~ add(add.sd) + prop(prop.sd) })
  }
  rC <- .rpemFit(rxode2::rxode2(modC), datC, ctlC)
  rR <- .rpemFit(rxode2::rxode2(modC), datC, ctlR)
  expect_equal(rC$addSd, rR$addSd, tolerance = 0.1)
  expect_equal(rC$propSd, rR$propSd, tolerance = 0.1)
  expect_equal(rC$addSd, 0.05, tolerance = 0.3)
  expect_equal(rC$propSd, 0.10, tolerance = 0.2)

  # TBS (boxCox): lambda is golden-sectioned in C++
  bcI <- function(y, l) if (l == 0) exp(y) else (l * y + 1)^(1 / l)
  datT <- mkdat(7, function(cp) {
    tcp <- ((pmax(cp, 1e-3))^0.5 - 1) / 0.5
    pmax(bcI(tcp + rnorm(length(cp), 0, 0.1), 0.5), 1e-3)
  })
  modT <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.2; lambda <- 1; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt()
            cp ~ add(add.sd) + boxCox(lambda) })
  }
  tC <- .rpemFit(rxode2::rxode2(modT), datT, ctlC)
  tR <- .rpemFit(rxode2::rxode2(modT), datT, ctlR)
  expect_equal(tC$addSd, tR$addSd, tolerance = 0.1)
  expect_equal(tC$lambda, tR$lambda, tolerance = 0.15)
})

test_that("est=rpem cLoop estimates a non-mu-ref structural beta (matches R loop)", {
  skip_on_cran()
  skip_on_ci()  # heavy: structural re-solve M-step in C++

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(9); nsub <- 40L; obsT <- seq(0.5, 24, by = 1.5); et0 <- rnorm(nsub, 0, sqrt(0.3))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = et0[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; d$DV <- 0; o <- d$evid == 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, 0.1)
    d
  }))
  # tcl is structural (no eta, not fixed); tv is fixed
  mod <- function() {
    ini({ tka <- 0.3; tcl <- 0.5; tv <- fix(3.45); add.sd <- 0.2; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  ui <- rxode2::rxode2(mod)
  ctlR <- rpemControl(nGauss = 200L, nMH = 60000L, mhBurn = 6000L, niter = 25L,
                      collect = 10L, seed = 42L, cores = 4L, cLoop = FALSE)
  ctlC <- rpemControl(nGauss = 200L, nMH = 60000L, mhBurn = 6000L, niter = 25L,
                      collect = 10L, seed = 42L, cores = 4L, cLoop = TRUE)
  rC <- .rpemFit(ui, dat, ctlC)
  rR <- .rpemFit(ui, dat, ctlR)
  # the structural beta (tcl) recovers the truth (1.0) and matches the R loop
  expect_equal(unname(rC$struct["tcl"]), 1.0, tolerance = 0.06)
  expect_equal(unname(rC$struct["tcl"]), unname(rR$struct["tcl"]), tolerance = 0.05)
  # reproducible run-to-run
  rC2 <- .rpemFit(ui, dat, ctlC)
  expect_equal(rC$struct, rC2$struct)
  expect_equal(rC$mu, rC2$mu)
})

test_that("est=rpem cLoop fits a mixture (matches R loop, reproducible)", {
  skip_on_cran()
  skip_on_ci()  # heavy: several mixture EM fits

  sim <- rxode2::rxode2({ ka <- exp(tka + eka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(52); nsub <- 120L; obsT <- seq(0.5, 24, by = 2)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    k <- if (stats::runif(1) < 0.5) 1L else 2L
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(sim, params = c(tka = c(0.1, 1.3)[k], tcl = 1, tv = 3.45,
                                         eka = stats::rnorm(1, 0, sqrt(0.3))),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
    d$DV[o] <- s$cp + stats::rnorm(nrow(s), 0, 0.1)
    d
  }))
  mod <- function() {
    ini({ tka1 <- 0.1; tka2 <- 1.3; tcl <- fix(1.0); tv <- fix(3.45)
          p1 <- 0.5; add.sd <- 0.2; eta.ka ~ 0.3 })
    model({ ka <- mix(exp(tka1 + eta.ka), p1, exp(tka2 + eta.ka))
            cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  ui <- rxode2::rxode2(mod)
  ctlR <- rpemControl(nGauss = 400L, nMH = 80000L, mhBurn = 8000L, niter = 30L,
                      collect = 10L, seed = 7L, cores = 4L, cLoop = FALSE)
  ctlC <- rpemControl(nGauss = 400L, nMH = 80000L, mhBurn = 8000L, niter = 30L,
                      collect = 10L, seed = 7L, cores = 4L, cLoop = TRUE)

  # full mixture fit object builds via the C++ loop
  fitC <- suppressMessages(nlmixr2(mod, dat, est = "rpem", control = ctlC))
  expect_s3_class(fitC, "nlmixr2FitData")
  expect_false(is.null(fitC$mixNum))

  rC <- .rpemFit(ui, dat, ctlC)
  rR <- .rpemFit(ui, dat, ctlR)
  # components separated and close to the R loop
  expect_lt(rC$mix$muK[1], 0.5)
  expect_gt(rC$mix$muK[2], 1.0)
  expect_equal(unname(rC$mix$muK[2]), unname(rR$mix$muK[2]), tolerance = 0.1)
  expect_equal(rC$addSd, rR$addSd, tolerance = 0.08)
  # reproducible run-to-run
  rC2 <- .rpemFit(ui, dat, ctlC)
  expect_equal(rC$mix$muK, rC2$mix$muK)
  expect_equal(rC$mix$w, rC2$mix$w)
})

test_that("est=rpem cLoop fits a mixture with a combined residual (matches R loop)", {
  skip_on_cran()
  skip_on_ci()  # heavy: mixture EM with the combined optimizer

  sim <- rxode2::rxode2({ ka <- exp(tka + eka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(42); nsub <- 120L; obsT <- seq(0.5, 24, by = 2.5)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    k <- if (stats::runif(1) < 0.6) 1L else 2L
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(sim, params = c(tka = c(0, 1.4)[k], tcl = 1, tv = 3.45,
                                         eka = stats::rnorm(1, 0, sqrt(0.3))),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
    d$DV[o] <- s$cp + stats::rnorm(nrow(s), 0, sqrt(0.05^2 + (0.1 * s$cp)^2))
    d
  }))
  mod <- function() {
    ini({ tka1 <- 0.2; tka2 <- 1.2; tcl <- fix(1.0); tv <- fix(3.45)
          p1 <- 0.5; add.sd <- 0.1; prop.sd <- 0.15; eta.ka ~ 0.3 })
    model({ ka <- mix(exp(tka1 + eta.ka), p1, exp(tka2 + eta.ka))
            cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) + prop(prop.sd) })
  }
  ui <- rxode2::rxode2(mod)
  ctlR <- rpemControl(nGauss = 250L, nMH = 60000L, mhBurn = 6000L, niter = 20L,
                      collect = 8L, seed = 7L, cores = 4L, cLoop = FALSE)
  ctlC <- rpemControl(nGauss = 250L, nMH = 60000L, mhBurn = 6000L, niter = 20L,
                      collect = 8L, seed = 7L, cores = 4L, cLoop = TRUE)
  cR <- .rpemFit(ui, dat, ctlR)
  cC <- .rpemFit(ui, dat, ctlC)
  # both residual parameters match the R loop and recover (true add 0.05, prop 0.10)
  expect_equal(cC$addSd, cR$addSd, tolerance = 0.1)
  expect_equal(cC$propSd, cR$propSd, tolerance = 0.1)
  expect_true(cC$addSd > 0.02 && cC$addSd < 0.09)
  expect_true(cC$propSd > 0.06 && cC$propSd < 0.14)
  # components separated
  expect_lt(cC$mix$muK[1], 0.6)
  expect_gt(cC$mix$muK[2], 1.0)
})

# cLoop covariate regression M-step (design/rpem/12 M5): a non-time-varying mu2
# covariate coefficient is estimated by the C++ weighted-regression M-step (nEta==1),
# rather than falling back to the R loop.  Must match the R loop and be reproducible.

test_that("est=rpem cLoop estimates a mu2 covariate via the C++ regression M-step", {
  skip_on_cran()
  skip_on_ci()  # heavy: multiple RPEM fits

  simMod <- rxode2::rxode2({ ka <- exp(0.45 + 0.35 * NTV + eka); cl <- exp(1)
                             v <- exp(3.45); cp <- linCmt() })
  set.seed(7); nsub <- 60L; obsT <- seq(0.5, 24, by = 1.5)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ntv <- rnorm(1, 0, 1); eka <- rnorm(1, 0, sqrt(0.12)); n <- length(obsT)
    ev <- data.frame(id = i, time = c(0, obsT), evid = c(1, rep(0, n)), amt = c(100, rep(0, n)),
                     cmt = 1, NTV = ntv, eka = eka)
    s <- rxode2::rxSolve(simMod, ev, returnType = "data.frame", addDosing = FALSE)
    ev$DV <- 0; o <- ev$evid == 0; ev$DV[o] <- s$cp + rnorm(sum(o), 0, 0.1); ev$eka <- NULL; ev
  }))
  mod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); lv <- fix(3.45); b_ntv <- 0.1; add.sd <- 0.2; eta.ka ~ 0.3 })
    model({ ka <- exp(tka + b_ntv * NTV + eta.ka); cl <- exp(tcl); v <- exp(lv)
            cp <- linCmt(); cp ~ add(add.sd) })
  }
  ui <- rxode2::rxode2(mod)
  ctlR <- rpemControl(nGauss = 300L, nMH = 60000L, mhBurn = 6000L, niter = 25L,
                      collect = 10L, seed = 1L, cores = 4L, cLoop = FALSE)
  ctlC <- rpemControl(nGauss = 300L, nMH = 60000L, mhBurn = 6000L, niter = 25L,
                      collect = 10L, seed = 1L, cores = 4L, cLoop = TRUE)
  rxode2::rxSetSeed(42); rfR <- .rpemFit(ui, dat, ctlR)
  rfC <- .rpemFit(ui, dat, ctlC)

  # the C++ regression M-step recovers the covariate coefficient (true 0.35) and
  # matches the R loop closely
  expect_true("b_ntv" %in% names(rfC$covCoef))
  expect_equal(unname(rfC$covCoef["b_ntv"]), 0.35, tolerance = 0.15)
  expect_equal(unname(rfC$covCoef["b_ntv"]), unname(rfR$covCoef["b_ntv"]), tolerance = 0.05)
  expect_equal(rfC$addSd, rfR$addSd, tolerance = 0.06)

  # reproducible run-to-run for a fixed core count
  rfC2 <- .rpemFit(ui, dat, ctlC)
  expect_equal(rfC$covCoef, rfC2$covCoef)
  expect_equal(rfC$mu, rfC2$mu)
})

# Dynamic-iteration seed stability (imp.cpp even/odd threefry streams): extending niter
# reproduces the exact per-iteration prefix of a shorter run at the same seed, so a phase
# of estimation can be lengthened without changing the shared iterations.

test_that("est=rpem cLoop is dynamic-iteration stable (longer run shares the prefix)", {
  skip_on_cran()
  skip_on_ci()  # heavy: two C++ loops

  struct <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(8); nsub <- 30L; etas <- rnorm(nsub, 0, sqrt(0.3)); obsT <- seq(0.5, 24, by = 2)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eta = etas[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, 0.1); d
  }))
  mod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.2; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  ui <- rxode2::rxUiDecompress(rxode2::rxode2(mod))
  cl <- .rpemClassify(ui)
  .nm <- c(paste0("THETA[", seq_len(cl$nTheta), "]"), paste0("ETA[", seq_len(cl$nEta), "]"))
  e <- new.env(); e$predOnly <- ui$rpemRxModel$predOnly
  e$rxControl <- rxode2::rxControl(atol = 1e-8, rtol = 1e-8, cores = 2L)
  e$param <- stats::setNames(cl$base, .nm); e$data <- dat
  runN <- function(ni) rpemEMLoopK1(e, cl$base, cl$etaIdx, cl$muIdx, cl$addSdIdx, cl$errType,
    cl$mu0, diag(as.matrix(cl$omega0)), cl$addSd0, c(-1L, -1L, -1L), c(0, 0, 0),
    as.integer(cl$structIdx), as.numeric(cl$struct0), ni, 200L, 2L, 30000L, 3000L, 123L,
    matrix(0, 0, 0), integer(0), numeric(0), numeric(0), integer(0), 0L, 8L, 5L, 1e7, 0, 20L, 1.0)
  r20 <- runN(20L); r30 <- runN(30L)
  expect_equal(r20$muTrace[1:20, , drop = FALSE], r30$muTrace[1:20, , drop = FALSE])
  expect_equal(r20$omegaTrace[1:20, , drop = FALSE], r30$omegaTrace[1:20, , drop = FALSE])
  expect_equal(as.numeric(r20$sdTrace)[1:20], as.numeric(r30$sdTrace)[1:20])
  expect_equal(as.numeric(r20$lnL)[1:20], as.numeric(r30$lnL)[1:20])
})

# Mode-centered importance sampling (impInflate) runs in the C++ cLoop too: the eta draw
# is centered at each subject's EBE (kept in C++) and importance-weighted, matching the R
# loop.  impInflate == 0 is byte-identical to prior sampling (logRatio == 0).

test_that("est=rpem cLoop supports mode-centered IS (impInflate) and matches the R loop", {
  skip_on_cran()
  skip_on_ci()  # heavy: several multi-eta RPEM fits

  sim <- rxode2::rxode2({ ka <- exp(lka + eka); cl <- exp(lcl + ecl); v <- exp(lv + ev); cp <- linCmt() })
  set.seed(42); nsub <- 50L; obsT <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 16, 24)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    eka <- rnorm(1, 0, sqrt(0.3)); ecl <- rnorm(1, 0, sqrt(0.1)); ev <- rnorm(1, 0, sqrt(0.06))
    ev0 <- rxode2::et(amt = 100, cmt = "depot"); ev0 <- rxode2::et(ev0, obsT)
    s <- rxode2::rxSolve(sim, params = c(lka = 0.5, lcl = 1, lv = 3.45, eka = eka, ecl = ecl, ev = ev),
                         events = ev0, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev0); d$id <- i; o <- d$evid == 0; d$DV <- 0
    d$DV[o] <- s$cp + rnorm(sum(o), 0, 0.1); d
  }))
  mod <- function() {
    ini({ lka <- 0.5; lcl <- fix(1.0); lv <- fix(3.45); eta.ka ~ 0.3; eta.cl ~ 0.1; eta.v ~ 0.06; add.sd <- 0.15 })
    model({ ka <- exp(lka + eta.ka); cl <- exp(lcl + eta.cl); v <- exp(lv + eta.v)
            cp <- linCmt(); cp ~ add(add.sd) })
  }
  ui <- rxode2::rxUiDecompress(rxode2::rxode2(mod))
  ctl <- function(ci, cl) rpemControl(nGauss = 300L, nMH = 50000L, mhBurn = 5000L, niter = 20L,
                                      collect = 8L, seed = 1L, cores = 4L, impInflate = ci, cLoop = cl)
  # impInflate=0: the C++ loop matches the R loop (prior sampling, logRatio 0)
  r0R <- .rpemFit(ui, dat, ctl(0, FALSE)); r0C <- .rpemFit(ui, dat, ctl(0, TRUE))
  expect_equal(r0C$omega[1], r0R$omega[1], tolerance = 0.02)
  # impInflate=4: mode-centering lifts the under-covered omega, and the C++ loop matches R
  r4R <- .rpemFit(ui, dat, ctl(4, FALSE)); r4C <- .rpemFit(ui, dat, ctl(4, TRUE))
  expect_equal(r4C$omega[1], r4R$omega[1], tolerance = 0.02)
  expect_gt(r4C$omega[1], r0C$omega[1])            # mode-centering raises om.ka in C++
})
