test_that("impmapControl option sanity and IS defaults", {
  expect_error(impmapControl(), NA)
  .ctl <- impmapControl()
  expect_s3_class(.ctl, "impmapControl")
  # inherits the FOCEI/MAP option surface
  expect_false(is.null(.ctl$maxOuterIterations))
  # mu-referencing is forced on for the MAP proposal center
  expect_identical(.ctl$muModel, "lin")
  # importance-sampling / EM defaults
  expect_identical(.ctl$isample, 300L)
  expect_identical(.ctl$nIter, 100L)
  expect_identical(.ctl$mapIter, 1L)
  expect_identical(.ctl$gamma, 1.0)
  expect_identical(.ctl$iscaleMin, 0.1)
  expect_identical(.ctl$iscaleMax, 10.0)
  expect_identical(.ctl$iaccept, 0.4)
  expect_null(.ctl$ctol)
  expect_identical(.ctl$nConvWindow, 10L)
  expect_identical(.ctl$impSeed, 42L)

  # round-trips through do.call (used by getValidNlmixrCtl / .foceiFamilyControl)
  expect_error(do.call(impmapControl, .ctl), NA)
})

test_that("impmapControl overrides and focei passthrough", {
  .ctl <- impmapControl(isample = 50L, gamma = 2, impSeed = 7L,
                        maxOuterIterations = 3L)
  expect_identical(.ctl$isample, 50L)
  expect_identical(.ctl$gamma, 2.0)
  expect_identical(.ctl$impSeed, 7L)
  # forwarded to foceiControl via ...
  expect_identical(.ctl$maxOuterIterations, 3L)
})

test_that("down-conversion to foceiControl strips IS-only names", {
  .env <- new.env()
  .env$impmapControl <- impmapControl()
  .fc <- nlmixr2est:::.impmapControlToFoceiControl(.env, assign = FALSE)
  expect_s3_class(.fc, "foceiControl")
  # IS/EM-only names must not leak into the plain foceiControl
  for (.n in nlmixr2est:::.impmapIsControlNames) {
    expect_null(.fc[[.n]])
  }
  # MAP-relevant focei options are preserved
  expect_identical(.fc$muModel, "lin")
})

test_that("impmap dispatch is registered and discoverable", {
  expect_true("impmap" %in% nlmixr2AllEst())
  expect_true(is.function(getS3method("nlmixr2Est", "impmap")))
  # mu-hook activation gate is a control-dependent predicate (mufocei-style)
  expect_true(is.function(attr(nlmixr2Est.impmap, "mu")))
})

test_that("getValidNlmixrCtl.impmap yields a default impmapControl", {
  expect_s3_class(getValidNlmixrCtl.impmap(list(NULL)), "impmapControl")
})

test_that("M1: impmap MAP pass exposes per-subject mode and Hessian", {
  one.cmt <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  .dat <- nlmixr2data::theo_sd
  # One EM iteration is enough to exercise the MAP + per-subject Hessian.
  .imp <- suppressWarnings(
    nlmixr2(one.cmt, .dat, "impmap", impmapControl(print = 0L, nIter = 1L)))

  expect_true(inherits(.imp, "nlmixr2FitCore"))
  # The MAP pass stashes each subject's mode + eta Hessian; check the Hessian is
  # present, square, symmetric, and positive-definite for subject 1.
  .env <- .imp$env
  expect_true(is.matrix(.env$impEtaMode) &&
                nrow(.env$impEtaMode) == length(unique(.dat$ID)))
  .H <- .env$impEtaHess
  expect_true(is.list(.H) && length(.H) == length(unique(.dat$ID)))
  .H1 <- .H[[1]]
  expect_true(is.matrix(.H1) && all(dim(.H1) == c(2, 2)))
  expect_equal(.H1, t(.H1), tolerance = 1e-6)
  expect_true(all(eigen(.H1, symmetric = TRUE, only.values = TRUE)$values > 0))
})

test_that("M2: threefry proposal sampler matches N(mode, gamma*H^-1)", {
  one.cmt <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  .dat <- nlmixr2data::theo_sd
  .nsub <- length(unique(.dat$ID))
  .gamma <- 1.5
  rxode2::rxSetSeed(42)
  .f <- suppressWarnings(
    nlmixr2(one.cmt, .dat, "impmap",
            impmapControl(print = 0L, nIter = 1L, isample = 4000L, gamma = .gamma)))
  .e <- .f$env
  expect_identical(.e$impNsample, 4000L)
  expect_equal(.e$impGammaUsed, .gamma)
  .S <- .e$impSamples
  expect_true(is.list(.S) && length(.S) == .nsub)
  expect_true(all(dim(.S[[1]]) == c(4000L, 2L)))
  # subject 1: empirical mean ~ mode, empirical cov ~ gamma * H^-1
  .m1 <- as.numeric(.e$impEtaMode[1, ])
  expect_equal(colMeans(.S[[1]]), .m1, tolerance = 0.03)
  .covTarget <- .gamma * solve(.e$impEtaHess[[1]])
  expect_equal(unname(cov(.S[[1]])), unname(.covTarget), tolerance = 0.02)
})

test_that("M2: sampler is thread-count independent (D6)", {
  one.cmt <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  .dat <- nlmixr2data::theo_sd
  .thr0 <- rxode2::getRxThreads()
  on.exit(rxode2::setRxThreads(.thr0), add = TRUE)
  .run <- function(nthr) {
    rxode2::setRxThreads(nthr)
    rxode2::rxSetSeed(42)
    suppressWarnings(
      nlmixr2(one.cmt, .dat, "impmap",
              impmapControl(print = 0L, nIter = 1L, isample = 100L)))$env$impSamples
  }
  .s1 <- .run(1L)
  .s4 <- .run(4L)
  # same base seed + per-subject reseed => bit-identical regardless of thread count
  for (.i in seq_along(.s1)) {
    expect_identical(.s1[[.i]], .s4[[.i]])
  }
})

test_that("M3: importance weights recover the conditional mean/variance", {
  one.cmt <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  .dat <- nlmixr2data::theo_sd
  .gamma <- 2
  rxode2::rxSetSeed(42)
  .f <- suppressWarnings(
    nlmixr2(one.cmt, .dat, "impmap",
            impmapControl(print = 0L, nIter = 1L, isample = 6000L, gamma = .gamma)))
  .e <- .f$env
  # E-step outputs present and well-formed
  expect_true(is.numeric(.e$impObj) && is.finite(.e$impObj))
  expect_equal(nrow(.e$impCondMean), length(.e$impNeff))
  expect_true(all(is.finite(.e$impLi)))
  # effective sample size strictly between 1 and isample
  expect_true(all(.e$impNeff > 1 & .e$impNeff <= 6000))

  .k <- 1L
  .mode1 <- as.numeric(.e$impEtaMode[.k, ])
  .H1 <- .e$impEtaHess[[.k]]
  # conditional mean tracks the mode for a near-Gaussian posterior
  expect_equal(as.numeric(.e$impCondMean[.k, ]), .mode1, tolerance = 0.03)
  # KEY: the importance weights reweight samples drawn from the inflated
  # proposal (cov = gamma*H^-1) back to the TRUE posterior covariance ~ H^-1,
  # NOT the proposal covariance.  This is what validates the weights.
  .B1 <- .e$impCondVar[[.k]]
  .postCov <- solve(.H1)
  .propCov <- .gamma * .postCov
  expect_equal(unname(.B1), unname(.postCov), tolerance = 0.02)
  # and B is clearly closer to H^-1 than to the proposal covariance
  expect_lt(max(abs(.B1 - .postCov)), max(abs(.B1 - .propCov)))
})

test_that("M4: EM converges to FOCEI on the mu-referenced params and Omega", {
  # Non-mu parameters (tv, add.sd) are held fixed so this isolates the EM update
  # of the mu-referenced thetas and Omega (the non-mu FD updates are a later
  # module); impmap should then match FOCEI on tka/tcl and the Omega diagonal.
  mfix <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- fix(3.45)
      eta.ka ~ 0.6; eta.cl ~ 0.3
      add.sd <- fix(0.7)
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  .d <- nlmixr2data::theo_sd
  .ff <- suppressWarnings(nlmixr2(mfix, .d, "focei", foceiControl(print = 0L, covMethod = "")))
  rxode2::rxSetSeed(42)
  .fi <- suppressWarnings(nlmixr2(mfix, .d, "impmap",
                                  impmapControl(print = 0L, nIter = 40L, isample = 300L)))
  expect_true(inherits(.fi, "nlmixr2FitCore"))
  expect_equal(fixef(.fi)[c("tka", "tcl")], fixef(.ff)[c("tka", "tcl")], tolerance = 0.05)
  expect_equal(unname(diag(.fi$omega)), unname(diag(.ff$omega)), tolerance = 0.1)
})

test_that("M4: mu-referenced covariate (updateMuGroups) is estimated", {
  # cl.wt is a mu-referenced covariate effect -- handled by the covariate
  # regression update (updateMuGroups), which impmap must drive.  The estimate
  # should match FOCEI, and the fit should report a nonzero mu covariate group.
  mcov <- function() {
    ini({
      tka <- 0.45; tcl <- 1; cl.wt <- 0.75; tv <- fix(3.45)
      eta.ka ~ 0.6; eta.cl ~ 0.3
      add.sd <- fix(0.7)
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + cl.wt * log(WT / 70) + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  .d <- nlmixr2data::theo_sd
  .ff <- suppressWarnings(nlmixr2(mcov, .d, "focei", foceiControl(print = 0L, covMethod = "")))
  rxode2::rxSetSeed(42)
  .fi <- suppressWarnings(nlmixr2(mcov, .d, "impmap",
                                  impmapControl(print = 0L, nIter = 40L, isample = 300L)))
  expect_true(inherits(.fi, "nlmixr2FitCore"))
  # a mu covariate group was actually set up and driven
  expect_true(.fi$env$impMuGroupN >= 1L)
  # the covariate coefficient matches FOCEI
  expect_equal(unname(fixef(.fi)["cl.wt"]), unname(fixef(.ff)["cl.wt"]), tolerance = 0.03)
  expect_equal(fixef(.fi)[c("tka", "tcl")], fixef(.ff)[c("tka", "tcl")], tolerance = 0.05)
})
