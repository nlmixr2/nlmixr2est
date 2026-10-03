# The one way a covariance is installed on a fit (R/cov.R): the guard, the
# covList stash, the condition refresh, and the setCov() paths built on them.

.pdCov <- function(nm = c("tka", "tcl")) {
  matrix(c(0.04, 0.01, 0.01, 0.09), 2, dimnames = list(nm, nm))
}

.fakeFitEnv <- function(cov = .pdCov(), covMethod = "r,s") {
  .e <- new.env(parent = emptyenv())
  .e$cov <- cov
  .e$covMethod <- covMethod
  .e$objDf <- data.frame(OBJF = 1, "Condition#(Cov)" = 99, "Condition#(Cor)" = 99, check.names = FALSE)
  .e
}

test_that(".covGuard() names why a matrix is rejected", {
  expect_identical(.covGuard(NULL)$reason, "could not be computed")
  expect_identical(.covGuard(matrix(1, 2, 3))$reason, "could not be computed")
  expect_identical(.covGuard(matrix(numeric(0), 0, 0))$reason, "could not be computed")
  expect_identical(.covGuard(matrix(c(1, NA, NA, 1), 2))$reason, "is not finite")
  expect_identical(.covGuard(matrix(c(1, Inf, Inf, 1), 2))$reason, "is not finite")
  expect_identical(.covGuard(matrix(c(1, 0.5, 0, 1), 2))$reason, "is not symmetric")
  expect_identical(.covGuard(matrix(c(1, 2, 2, 1), 2))$reason, "is not positive definite")
  expect_identical(.covGuard(matrix(c(0, 0, 0, 1), 2))$reason, "is not positive definite")
  expect_identical(.covGuard(matrix(0, 2, 2))$reason, "is not positive definite")
  expect_false(.covGuard(matrix(c(1, 2, 2, 1), 2))$ok)
  # a finite square matrix comes back symmetrized even when rejected
  expect_equal(.covGuard(matrix(c(1, 2, 2, 1), 2))$cov, matrix(c(1, 2, 2, 1), 2))
  expect_null(.covGuard(matrix(c(1, NA, NA, 1), 2))$cov)
})

test_that(".covGuard() accepts a positive-definite matrix and symmetrizes it", {
  .m <- .pdCov()
  .m[1, 2] <- .m[1, 2] + 1e-12 # rounding-level asymmetry is accepted
  .g <- .covGuard(.m)
  expect_true(.g$ok)
  expect_identical(.g$reason, "")
  expect_identical(.g$cov[1, 2], .g$cov[2, 1])
  expect_equal(.g$cov[1, 2], 0.01 + 5e-13, tolerance = 1e-15)
  expect_identical(dimnames(.g$cov), dimnames(.m))
  expect_equal(.g$ev, eigen(.g$cov, symmetric = TRUE, only.values = TRUE)$values)
})

test_that(".covMethodFromSlot() inverts the foceiControl() covMethod slot", {
  expect_identical(.covMethodFromSlot(1L), "r,s")
  expect_identical(.covMethodFromSlot(2L), "r")
  expect_identical(.covMethodFromSlot(3L), "s")
  expect_identical(.covMethodFromSlot(2, "analytic"), "analytic")
  expect_identical(.covMethodFromSlot(1L, "analytic"), "r,s")
  expect_identical(.covMethodFromSlot(0L), "")
  expect_identical(.covMethodFromSlot(NA_integer_), "")
  expect_identical(.covMethodFromSlot("r"), "")
  expect_identical(.covMethodFromSlot(NULL), "")
  for (.n in names(.covMethodSlot)) {
    expect_identical(.covMethodFromSlot(.covMethodSlot[[.n]]), .n)
  }
})

test_that("a control given a foceiControl() covMethod slot keeps the covariance the slot names", {
  for (.n in names(.covMethodSlot)) {
    .slot <- foceiControl(covMethod = .n)$covMethod
    expect_identical(nlmeControl(covMethod = .slot)$covMethod, .n)
    expect_identical(vaeControl(covMethod = .slot)$covMethod, .n)
    expect_identical(emviControl(covMethod = .slot)$covMethod, .n)
  }
  expect_identical(nlmeControl(covMethod = 0L)$covMethod, "")
  expect_identical(vaeControl(covMethod = 0)$covMethod, "")
  expect_error(nlmeControl(covMethod = 4L), "foceiControl() slot", fixed = TRUE)
})

test_that("a covariance refit sets each option under its own name only", {
  .obj <- new.env(parent = emptyenv())
  .obj$foceiControl <- foceiControl()
  local_mocked_bindings(
    getData = function(object) NULL,
    nlmixr2CreateOutputFromUi = function(ui, data, control, ...) control
  )
  # hessEps, rmatNorm and gillStepCov are prefixes of the log-likelihood options
  .ctl <- .setCovRefit(.obj, hessEpsLlik = 1e-3, rmatNormLlik = 0L, gillStepCovLlik = 3)
  expect_identical(
    .ctl[c("hessEpsLlik", "rmatNormLlik", "gillStepCovLlik")],
    list(hessEpsLlik = 1e-3, rmatNormLlik = 0L, gillStepCovLlik = 3)
  )
  expect_identical(
    .ctl[c("hessEps", "rmatNorm", "gillStepCov")],
    .obj$foceiControl[c("hessEps", "rmatNorm", "gillStepCov")]
  )
})

test_that(".covInnerIterations() gives a refit's covariance legs the fit's inner budget", {
  expect_identical(.covInnerIterations(foceiControl(maxInnerIterations = 250L)), 250L)
  # a fit with no budget of its own (saem, nlme, vae and vi evaluate their ETAs) gets
  # foceiControl()'s default
  .def <- 1000L
  expect_identical(as.integer(formals(foceiControl)$maxInnerIterations), .def)
  expect_identical(.covInnerIterations(foceiControl(maxInnerIterations = 0L)), .def)
  expect_identical(.covInnerIterations(list()), .def)
  expect_identical(.covInnerIterations(list(maxInnerIterations = NA_integer_)), .def)
  expect_identical(.covInnerIterations(list(maxInnerIterations = -3L)), .def)
  expect_identical(.covInnerIterations(list(maxInnerIterations = Inf)), .def)
  expect_identical(.covInnerIterations(list(maxInnerIterations = c(5L, 6L))), .def)
})

test_that("a covariance refit reports the fit's ETAs but differentiates its marginal likelihood", {
  .obj <- new.env(parent = emptyenv())
  .obj$foceiControl <- foceiControl(maxInnerIterations = 250L, interaction = TRUE)
  local_mocked_bindings(
    getData = function(object) NULL,
    nlmixr2CreateOutputFromUi = function(ui, data, control, ...) control
  )
  .ctl <- .setCovRefit(.obj, covMethod = "r")
  # the refit evaluates the fit's ETAs; its covariance legs optimize them with the fit's
  # budget, and with the fit's interaction
  expect_identical(.ctl$maxInnerIterations, 0L)
  expect_identical(.ctl$covMaxInnerIterations, 250L)
  expect_identical(.ctl$interaction, 1L)
})

test_that("covMaxInnerIterations survives a foceiControl() round trip", {
  # nlmixr2() rebuilds a control it is given (getValidNlmixrCtl()), and the post-fit
  # recompute hands its control to nlmixr2()
  .ctl <- foceiControl()
  expect_null(.ctl$covMaxInnerIterations)
  .ctl$covMaxInnerIterations <- 300L
  expect_identical(getValidNlmixrCtl.focei(list(.ctl))$covMaxInnerIterations, 300L)
  expect_identical(foceiControl(covMaxInnerIterations = 7)$covMaxInnerIterations, 7L)
  expect_error(foceiControl(covMaxInnerIterations = 0L), "covMaxInnerIterations")
  expect_error(foceiControl(covMaxInnerIterations = 1.5), "covMaxInnerIterations")
})

test_that(".covInstall() installs, stashes the replaced covariance and refreshes its diagnostics", {
  .e <- .fakeFitEnv()
  .new <- .pdCov() * 4
  expect_true(.covInstall(.e, .new, "sa"))
  expect_identical(.e$covMethod, "sa")
  expect_equal(.e$cov, .new)
  expect_identical(names(.e$covList), "r,s")
  expect_equal(.e$covList[["r,s"]], .pdCov())
  .ev <- eigen(.new, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(.e$objDf[["Condition#(Cov)"]], max(.ev) / min(.ev))
  expect_equal(.e$conditionNumberCov, max(.ev) / min(.ev))
  expect_equal(.e$eigenCov, .ev)
  expect_equal(unname(diag(.e$fullCor)), unname(sqrt(diag(.new))))
  # an unlabelled covariance is kept as "prev"; the installed name is never stashed
  .e2 <- .fakeFitEnv(covMethod = "")
  .covInstall(.e2, .new, "sa")
  expect_identical(names(.e2$covList), "prev")
  .covInstall(.e2, .new * 2, "sa")
  expect_identical(names(.e2$covList), "prev")
})

test_that(".covInstall() leaves the fit as it was when the matrix is rejected", {
  .e <- .fakeFitEnv()
  .bad <- matrix(c(1, 2, 2, 1), 2, dimnames = dimnames(.pdCov()))
  expect_warning(
    .ok <- .covInstall(.e, .bad, "sa", extras = list(parFixedDf = "refit")),
    "\"sa\" covariance is not positive definite; kept \"r,s\"",
    fixed = TRUE
  )
  expect_false(.ok)
  expect_identical(attr(.ok, "reason"), "is not positive definite")
  expect_identical(.e$covMethod, "r,s")
  expect_equal(.e$cov, .pdCov())
  expect_false(exists("covList", envir = .e, inherits = FALSE))
  expect_false(exists("parFixedDf", envir = .e, inherits = FALSE))
  expect_identical(.e$objDf[["Condition#(Cov)"]], 99)
  .e2 <- new.env(parent = emptyenv())
  expect_warning(
    .covInstall(.e2, NULL, "imp"),
    "\"imp\" covariance could not be computed; none installed",
    fixed = TRUE
  )
  expect_silent(.covInstall(.e2, NULL, "imp", warn = FALSE))
  expect_warning(.covInstall(.e2, NULL, NULL), "the covariance could not be computed; none installed", fixed = TRUE)
})

test_that(".covInstall() says when what it installed is not what was asked for", {
  .e <- .fakeFitEnv()
  expect_warning(
    .covInstall(.e, .pdCov(), "linFim", what = "sa"),
    "\"linFim\" covariance installed instead of the requested \"sa\"",
    fixed = TRUE
  )
  expect_identical(.e$covMethod, "linFim")
  # a correction decoration or the full-shape suffix names the same covariance
  expect_silent(.covInstall(.e, .pdCov(), "|r|,s (full)", what = "r,s"))
  expect_silent(.covInstall(.e, .pdCov(), "analytic (full)", what = "analytic"))
})

test_that("the post-fit mu recompute keeps the fit's covariance when its result is unusable", {
  .e <- .fakeFitEnv(covMethod = "nlme")
  local_mocked_bindings(.foceiRecomputeMuCov = function(fit, est) {
    list(
      cov = matrix(c(1, 2, 2, 1), 2, dimnames = dimnames(.pdCov())),
      covMethod = "r,s",
      extras = list(parFixedDf = "refit"),
      what = "r,s"
    )
  })
  expect_warning(
    .foceiInstallMuCov(.e, "nlme"),
    "\"r,s\" covariance is not positive definite; kept \"nlme\"",
    fixed = TRUE
  )
  expect_identical(.e$covMethod, "nlme")
  expect_equal(.e$cov, .pdCov())
  expect_false(exists("parFixedDf", envir = .e, inherits = FALSE))
  # a recompute that failed outright says so too
  local_mocked_bindings(.foceiRecomputeMuCov = function(fit, est) list(what = "analytic"))
  expect_warning(
    .foceiInstallMuCov(.e, "nlme"),
    "\"analytic\" covariance could not be computed; kept \"nlme\"",
    fixed = TRUE
  )
  expect_identical(.e$covMethod, "nlme")
  expect_equal(.e$cov, .pdCov())
  # no recompute requested: nothing to say
  local_mocked_bindings(.foceiRecomputeMuCov = function(fit, est) NULL)
  expect_silent(.foceiInstallMuCov(.e, "nlme"))
})

test_that("the post-fit mu recompute carries the refit's tables and stashes the fit's covariance", {
  .e <- .fakeFitEnv(covMethod = "nlme")
  local_mocked_bindings(.foceiRecomputeMuCov = function(fit, est) {
    list(
      cov = .pdCov() * 2,
      covMethod = "r,s (full)",
      extras = list(parFixedDf = "refit", covList = list(r = .pdCov())),
      what = "r,s"
    )
  })
  expect_true(.foceiInstallMuCov(.e, "nlme"))
  expect_identical(.e$covMethod, "r,s (full)")
  expect_equal(.e$cov, .pdCov() * 2)
  expect_identical(.e$parFixedDf, "refit")
  expect_identical(names(.e$covList), c("r", "nlme"))
  expect_equal(.e$covList$nlme, .pdCov())
})

nmTest({
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

  .fitOnce <- function() {
    suppressWarnings(nlmixr2(
      .oneCmt,
      nlmixr2data::theo_sd,
      est = "focei",
      control = foceiControl(print = 0, covMethod = "r,s", calcTables = FALSE)
    ))
  }

  test_that("a cached covariance that is not positive definite is not reinstalled", {
    .fit <- .fitOnce()
    .cov0 <- .fit$cov
    .m0 <- .fit$covMethod
    .bad <- .cov0
    .bad[1, 1] <- -.bad[1, 1]
    .cl <- .fit$env$covList
    .cl$imp <- .bad
    assign("covList", .cl, envir = .fit$env)
    expect_error(
      suppressMessages(setCov(.fit, "imp")),
      "covMethod=\"imp\" could not be computed for this fit (the result is not positive definite)",
      fixed = TRUE
    )
    expect_identical(.fit$covMethod, .m0)
    expect_equal(.fit$cov, .cov0)
  })

  test_that("setCov() checks the covariance a nested fit gives before installing it", {
    .fit <- .fitOnce()
    .cov0 <- .fit$cov
    .m0 <- .fit$covMethod
    local_mocked_bindings(.covRecomputeNative = function(fit, est, control, useEtaMat = TRUE) {
      list(cov = get("cov", envir = fit$env) * 2, covMethod = "linFim", mixRotated = TRUE)
    })
    expect_error(
      setCov(.fit, "sa"),
      "covMethod=\"sa\" could not be computed for this fit (the result was \"linFim\")",
      fixed = TRUE
    )
    expect_identical(.fit$covMethod, .m0)
    expect_equal(.fit$cov, .cov0)
    expect_null(.fit$env$covOptions$sa)
  })

  test_that("a finite-difference setCov() keeps the label its refit computed", {
    .fit <- .fitOnce()
    .theta <- .fit$env$covList[["r,s"]]
    .pf <- .fit$parFixedDf
    local_mocked_bindings(nlmixr2CreateOutputFromUi = function(...) {
      list(cov = .theta, covMethod = "|r|,s", parFixedDf = .pf, parFixed = NULL)
    })
    suppressMessages(setCov(.fit, "r,s", control = rsControl(hessEps = 1e-4)))
    expect_identical(.fit$covMethod, "|r|,s")
    expect_equal(.fit$cov, .theta)
    # the options are recorded under the installed name too, so the same request is not recomputed
    expect_equal(.fit$env$covOptions[["|r|,s"]]$hessEps, 1e-4)
    expect_error(setCov(.fit, "r,s", control = rsControl(hessEps = 1e-4)), "no need to switch")
  })

  test_that("a finite-difference setCov() whose refit falls back to another covariance installs nothing", {
    .fit <- .fitOnce()
    .cov0 <- .fit$cov
    .m0 <- .fit$covMethod
    .theta <- .fit$env$covList[["r,s"]]
    .pf <- .fit$parFixedDf
    local_mocked_bindings(nlmixr2CreateOutputFromUi = function(...) {
      list(cov = .theta, covMethod = "r", parFixedDf = .pf, parFixed = NULL)
    })
    expect_error(
      setCov(.fit, "r,s", control = rsControl(hessEps = 1e-4)),
      "covMethod=\"r,s\" could not be computed for this fit (the result was \"r\")",
      fixed = TRUE
    )
    expect_identical(.fit$covMethod, .m0)
    expect_equal(.fit$cov, .cov0)
  })

  test_that("a refit with no covariance (a mu-referenced fit's covariance step) installs nothing", {
    .fit <- .fitOnce()
    .cov0 <- .fit$cov
    .m0 <- .fit$covMethod
    local_mocked_bindings(nlmixr2CreateOutputFromUi = function(...) {
      list(cov = NULL, covMethod = "", parFixedDf = NULL, parFixed = NULL)
    })
    expect_error(
      setCov(.fit, "s", control = rsControl(hessEps = 1e-4)),
      "covMethod=\"s\" could not be computed for this fit; the covariance is left unchanged",
      fixed = TRUE
    )
    expect_identical(.fit$covMethod, .m0)
    expect_equal(.fit$cov, .cov0)
    # getVarCov(force = TRUE) refits the same way, and warns instead
    expect_warning(
      .v <- nlme::getVarCov(.fit, force = TRUE),
      sprintf("the covariance could not be computed; kept \"%s\"", .m0),
      fixed = TRUE
    )
    expect_equal(.v, .cov0)
    expect_identical(.fit$covMethod, .m0)
  })

  # The routes that compute a FOCEi covariance at a finished fit's estimates --
  # setCov(), getVarCov(force = TRUE), the post-fit recompute and the vae/vi
  # hand-off -- report the fit's ETAs but differentiate its marginal likelihood: they
  # give the estimation-time covariance of a FOCEi fit with no outer iterations
  # started from the fit's ETAs.  The full nlmixr2() refits (the recompute, vae/vi)
  # reproduce it exactly; setCov()'s lighter refit to within the inner problem's
  # tolerance (theta SEs within 1e-5 relative, Omega SEs within 2e-3 of the full
  # shape).  They used to hold the ETAs fixed (and setCov() dropped the
  # interaction), which made the theta SEs several times too small.
  .zeroOuterRef <- function(fit, covMethod, covFull = FALSE, interaction = TRUE) {
    suppressWarnings(nlmixr2(
      fit$finalUi,
      nlmixr2data::theo_sd,
      est = "focei",
      control = foceiControl(
        print = 0,
        calcTables = FALSE,
        maxOuterIterations = 0L,
        maxInnerIterations = 1000L,
        etaMat = fit$etaMat,
        covMethod = covMethod,
        covFull = covFull,
        interaction = interaction
      )
    ))
  }
  .seOf <- function(fit) sqrt(diag(fit$cov))
  .maxRel <- function(a, b) max(abs(a / b[names(a)] - 1))
  .ceOneCmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1.0
      tv <- 3.45
      add.sd <- 0.3
      prop.sd <- 0.1
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd) + prop(prop.sd)
    })
  }

  test_that("a finite-difference setCov() is the covariance of the fit's marginal likelihood", {
    .none <- suppressWarnings(nlmixr2(
      .oneCmt,
      nlmixr2data::theo_sd,
      est = "focei",
      control = foceiControl(print = 0, covMethod = "", calcTables = FALSE)
    ))
    .th <- c("tka", "tcl", "tv", "add.sd")
    for (.m in c("r", "s", "r,s", "r (full)")) {
      .full <- .covIsFull(.m)
      .ref <- .zeroOuterRef(.none, .covBaseName(.m), .full)
      expect_identical(.ref$covMethod, .m)
      suppressMessages(suppressWarnings(setCov(.none, .m)))
      expect_identical(.none$covMethod, .m)
      expect_setequal(names(.seOf(.none)), names(.seOf(.ref)))
      expect_lt(.maxRel(.seOf(.none)[.th], .seOf(.ref)), if (.full) 1e-3 else 1e-4, label = .m)
      if (.full) {
        .om <- setdiff(names(.seOf(.ref)), .th)
        expect_lt(.maxRel(.seOf(.none)[.om], .seOf(.ref)), 1e-2, label = .m)
      }
    }
    # getVarCov(force = TRUE) recomputes the same way
    .fit <- suppressWarnings(nlmixr2(
      .oneCmt,
      nlmixr2data::theo_sd,
      est = "focei",
      control = foceiControl(print = 0, covMethod = "r", covFull = FALSE, calcTables = FALSE)
    ))
    .v <- suppressMessages(suppressWarnings(nlme::getVarCov(.fit, force = TRUE)))
    expect_lt(.maxRel(sqrt(diag(.v)), .seOf(.zeroOuterRef(.fit, "r"))), 1e-4)
  })

  test_that("setCov() differentiates the likelihood the fit used: interaction is kept", {
    # with a proportional error the interaction changes the objective; setCov() used to
    # differentiate the FOCE one (interaction = 0) of a FOCEI fit
    .none <- suppressWarnings(nlmixr2(
      .ceOneCmt,
      nlmixr2data::theo_sd,
      est = "focei",
      control = foceiControl(print = 0, covMethod = "", calcTables = FALSE)
    ))
    .focei <- .zeroOuterRef(.none, "r")
    .foce <- .zeroOuterRef(.none, "r", interaction = FALSE)
    # the two likelihoods give tka SEs 30% apart here
    expect_gt(abs(.seOf(.foce)[["tka"]] / .seOf(.focei)[["tka"]] - 1), 0.2)
    suppressMessages(suppressWarnings(setCov(.none, "r")))
    expect_identical(.none$covMethod, "r")
    expect_lt(.maxRel(.seOf(.none), .seOf(.focei)), 1e-4)
  })

  test_that("the post-fit mu recompute differentiates the marginal likelihood", {
    .ctl <- foceiControl(print = 0, calcTables = FALSE)
    .mf <- suppressWarnings(nlmixr2(.oneCmt, nlmixr2data::theo_sd, est = "mfocei", control = .ctl))
    expect_identical(.mf$covMethod, "r,s (full)")
    .ref <- .zeroOuterRef(.mf, "r,s", covFull = TRUE)
    expect_setequal(names(.seOf(.mf)), names(.seOf(.ref)))
    expect_lt(.maxRel(.seOf(.mf), .seOf(.ref)), 1e-6)
    # an "analytic" request outside the analytic scope (linCmt) falls back to the
    # finite-difference sandwich, which is marginal too, and says so
    .an <- foceiControl(print = 0, calcTables = FALSE, covMethod = "analytic")
    expect_warning(
      .mfa <- nlmixr2(.oneCmt, nlmixr2data::theo_sd, est = "mfocei", control = .an),
      "\"r,s (full)\" covariance installed instead of the requested \"analytic\"",
      fixed = TRUE
    )
    expect_identical(.mfa$covMethod, "r,s (full)")
    expect_lt(.maxRel(.seOf(.mfa), .seOf(.ref)), 1e-6)
  })

  test_that("vae and emvi report the covariance of their likelihood's marginal at their estimates", {
    skip_on_cran()
    # Their output step evaluates the FOCE objective at their own ETAs; the covariance
    # comes after it, with the interaction of their `likelihood` (FOCEI) and the ETAs
    # optimized from theirs at every finite-difference leg.  It held the encoder or
    # variational means fixed, with interaction = 0.
    .refOf <- function(fit) {
      .ctl <- fit$foceiControl
      .ctl$maxInnerIterations <- 1000L
      .ctl$etaMat <- fit$etaMat
      .ctl$calcTables <- FALSE
      suppressMessages(suppressWarnings(nlmixr2(fit$finalUi, nlmixr2data::theo_sd, est = "focei", control = .ctl)))
    }
    .fits <- list(
      vae = vaeControl(
        itersBurnIn = 8L,
        iters = 16L,
        klWarmup = 4L,
        gammaIter = 12L,
        covariateSelection = FALSE,
        print = 0L,
        covMethod = "r",
        seed = 1L
      ),
      emvi = emviControl(iters = 60L, print = 0L, covMethod = "r", seed = 7L)
    )
    for (.est in names(.fits)) {
      .fit <- suppressMessages(suppressWarnings(nlmixr2(
        .oneCmt,
        nlmixr2data::theo_sd,
        est = .est,
        control = .fits[[.est]]
      )))
      # the fit's FOCEi control carries its likelihood's interaction; the output step's
      # control (in finalUi) evaluated FOCE
      expect_identical(.fit$foceiControl$interaction, 1L, label = .est)
      .ref <- .refOf(.fit)
      expect_identical(.covFdType(.fit$covMethod), "r", label = .est)
      expect_identical(.fit$covMethod, .ref$covMethod, label = .est)
      expect_setequal(names(.seOf(.fit)), names(.seOf(.ref)))
      expect_lt(.maxRel(.seOf(.fit), .seOf(.ref)), 1e-6, label = .est)
      # the method's own parameter table is kept, its SEs taken from the installed matrix
      .pf <- .fit$parFixedDf
      .th <- c("tka", "tcl", "tv", "add.sd")
      expect_setequal(rownames(.pf), .th)
      expect_equal(.pf[.th, "SE"], unname(.seOf(.fit)[.th]), tolerance = 1e-12, label = .est)
    }
  })
})
