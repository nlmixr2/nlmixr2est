nmTest({
  test_that("cov-focei", {
    dat <-
      warfarin |>
      dplyr::filter(dvid == "cp")

    #### doesn't work with FOCEI ### doesn't work with SAEM with CRWES=True but does iwth CWRES=FALSE
    One.SD.ODE <- function() {
      ini({
        # Where initial conditions/variables are specified
        lcl <- log(0.135) # log Cl (L/h)
        lv <- log(8) # log V (L)
        lmtt <- log(1.1) # log MTx

        prop.err <- 0.15 # proportional error (SD/mean)
        add.err <- 0.6 # additive error (mg/L)
        eta.cl + eta.v + eta.mtt ~ c(
          0.1,
          0.001, 0.1,
          0.001, 0.001, 0.1
        )
      })
      model({
        # Where the model is specified
        cl <- exp(lcl + eta.cl)
        v <- exp(lv + eta.v)
        mtt <- exp(lmtt + eta.mtt)
        ktr <- 6 / mtt

        ## ODE example
        d/dt(depot) <- -ktr * depot
        d/dt(central) <- ktr * trans5 - (cl / v) * central
        d/dt(trans1) <- ktr * depot - ktr * trans1
        d/dt(trans2) <- ktr * trans1 - ktr * trans2
        d/dt(trans3) <- ktr * trans2 - ktr * trans3
        d/dt(trans4) <- ktr * trans3 - ktr * trans4
        d/dt(trans5) <- ktr * trans4 - ktr * trans5

        ## Concentration is calculated
        cp <- central / v
        ## And is assumed to follow proportional and additive error
        cp ~ prop(prop.err) + add(add.err)
      })
    }

    f <- .nlmixr(One.SD.ODE, dat, "posthoc")

    expect_s3_class(f, "nlmixr2.posthoc")
  })

  test_that("covMethod r and s SEs are not inflated vs the sandwich (issue #666)", {
    one.cmt <- function() {
      ini({
        tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    # covFull = FALSE pins this to the NATIVE theta-only covariance, which is the
    # thing the ratios below are about.  With the FD-full covariance installed the
    # theta SEs come from inverting the JOINT theta+omega matrix instead, a
    # different quantity (it carries the omega estimation uncertainty), and
    # se_s/se_rs runs ~5 -- into the range this test reads as the old constant
    # factor returning.  Keep the guard on the estimator it was written for; the
    # SE/cov invariant under covFull is asserted separately below.
    fit_r <- .nlmixr(one.cmt, theo_sd, "focei", foceiControl(covMethod = "r", print = 0, covFull = FALSE))
    fit_s <- .nlmixr(one.cmt, theo_sd, "focei", foceiControl(covMethod = "s", print = 0, covFull = FALSE))
    fit_rs <- .nlmixr(one.cmt, theo_sd, "focei", foceiControl(covMethod = "r,s", print = 0, covFull = FALSE))

    p <- c("tka", "tcl", "tv")
    se_r <- fit_r$parFixedDf[p, "SE"]
    se_s <- fit_s$parFixedDf[p, "SE"]
    se_rs <- fit_rs$parFixedDf[p, "SE"]

    # Before the #666 fix: SE_r was ~sqrt(2)*SE_rs (covR was 2*Rinv) and SE_s was ~2x that again
    # (covS was 4*Sinv).  covMethod="r" (observed information) and the "r,s" sandwich obey the
    # information equality, so r/rs stays ~1 and a return of the sqrt(2) scaling would push it >1.41.
    expect_true(all(se_r / se_rs < 1.3), label = "covMethod='r' SE not inflated vs sandwich")
    # covMethod="s" is the OPG (score cross-product) estimator.  On this 12-subject dataset it
    # legitimately disagrees with the Hessian in its off-diagonals (S corr(tcl,tv) ~0.85 vs R's
    # ~0.12), so s/rs runs ~2 for cl/v -- finite-sample OPG behaviour that collapses toward 1 on
    # larger / well-specified data, NOT a scaling error.  This leg guards against the old constant
    # factor (covS = 4*Sinv), which would DOUBLE these ratios to ~4.5.
    #
    # Bound recalibrated from 3 to 3.5 when trustFterm/trustMterm's default was
    # tightened to 10^(-sigdig-2): tcl's ratio moved from 1.00 to 3.09.  The
    # looser inner solve was the reason the old numbers looked tame -- it also
    # put se_r/se_rs at 0.30-0.50, which CONTRADICTS the information equality
    # the r leg above asserts; tightening moves that to 0.88-1.26, i.e. onto the
    # ~1 this test says to expect.  3.5 still leaves clear room under the ~4.5 a
    # return of the constant factor would produce, which is what this guards.
    expect_true(all(se_s / se_rs < 3.5), label = "covMethod='s' SE not inflated by the old 2x constant factor")
  })

  test_that("reported SE matches sqrt(diag(fit$cov)) (nlmixr2extra#125)", {
    # .foceiInstallFdFullCov() replaces $cov after the C++ step has already
    # derived popDf$SE from the covariance it discards.  Without a refresh the
    # fit reports SEs describing a matrix it no longer holds, and a setCov()
    # round trip silently changes them.
    one.cmt <- function() {
      ini({
        tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    for (.full in c(TRUE, FALSE)) {
      .fit <- .nlmixr(one.cmt, theo_sd, "focei", foceiControl(print = 0, covFull = .full))
      .p <- intersect(rownames(.fit$parFixedDf), rownames(.fit$cov))
      expect_true(length(.p) > 0)
      expect_equal(
        unname(.fit$parFixedDf[.p, "SE"]),
        unname(sqrt(diag(.fit$cov))[.p]),
        label = paste0("covFull=", .full, " SE == sqrt(diag(cov))")
      )
    }
  })

  test_that("covariance with many omegas fixed will not crash focei", {
    one.compartment <- function() {
      ini({
        tka <- fix(0.45) # Log Ka
        tcl <- 1 # Log Cl
        tv <- fix(3.45)    # Log V
        eta.ka ~ fix(0.6)
        eta.cl ~ fix(0.3)
        eta.v ~ fix(0.1)
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.sd)
      })
    }
    fit <- .nlmixr(one.compartment, theo_sd, est = "focei", control = foceiControl(print = 0, maxOuterIterations = 0L))
    expect_s3_class(fit, "nlmixr2FitCore")
  })

  test_that("shi21maxOuter runs no step search of its own in the covariance step", {
    one.cmt <- function() {
      ini({
        tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    # innerOpt = "trust" counts every inner solve, so nTrustInner counts the
    # objective evaluations of the whole fit
    .ctl <- foceiControl(
      print = 0,
      maxOuterIterations = 0L,
      covFull = FALSE,
      innerOpt = "trust"
    )
    .gill <- .nlmixr(one.cmt, theo_sd, "focei", .ctl)
    .ctl$shi21maxOuter <- 8L
    .shi <- .nlmixr(one.cmt, theo_sd, "focei", .ctl)
    # The covariance steps are Gill's either way.  A Shi21 search used to run
    # first, be overwritten, and leave the inner problem where its last probe
    # put it.
    expect_identical(.shi$nTrustInner, .gill$nTrustInner)
    expect_identical(.shi$scaleInfo, .gill$scaleInfo)
    expect_identical(.shi$cov, .gill$cov)
  })

  test_that("every covariance stage is taken about the estimates", {
    one.cmt <- function() {
      ini({
        tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    # "s" runs S and the full FD straight after the step search.  "r,s" runs
    # the R stencil first, whose last leg used to leave theta at
    # theta0 - 2*eps0, where S then took its centre, and whose own last leg
    # moved it again before the full FD read it.
    .rs <- .nlmixr(one.cmt, theo_sd, "focei", foceiControl(print = 0))
    .s <- .nlmixr(one.cmt, theo_sd, "focei", foceiControl(print = 0, covMethod = "s"))
    expect_equal(.rs$env$S0, .s$env$S0, tolerance = 1e-6)
    expect_equal(.rs$env$.fdFullCov, .s$env$.fdFullCov, tolerance = 1e-6)
    expect_equal(.rs$env$.fdFullS, .s$env$.fdFullS, tolerance = 1e-6)
  })

  test_that("llikObs is that of the estimates, not of the last covariance leg", {
    one.cmt <- function() {
      ini({
        tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .ctl <- list(focei = foceiControl, laplace = laplaceControl)
    for (.est in names(.ctl)) {
      .def <- .nlmixr(one.cmt, theo_sd, .est, .ctl[[.est]](print = 0, maxOuterIterations = 0L))
      .none <- .nlmixr(
        one.cmt,
        theo_sd,
        .est,
        .ctl[[.est]](print = 0, maxOuterIterations = 0L, covMethod = "")
      )
      expect_false(is.null(.def$cov))
      expect_identical(.def$llikObs, .none$llikObs, label = .est)
    }
    # with no etas the objective is -2 * the sum of llikObs, and the covariance
    # step is the R matrix alone
    noEta <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka); cl <- exp(tcl); v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .pop <- .nlmixr(noEta, theo_sd, "focei", foceiControl(print = 0, maxOuterIterations = 0L))
    expect_false(is.null(.pop$cov))
    expect_equal(-2 * sum(.pop$llikObs, na.rm = TRUE), .pop$objf)
  })

  test_that("an agq fit reports llikObs at its ETAs, not at the last quadrature node", {
    # Every quadrature node re-solves the subject at another eta, rewriting its
    # per-observation log-likelihoods; the fit has to report them at the mode,
    # where the table's IPRED is solved.  No covariance step, so nothing else
    # moves them.
    one.cmt <- function() {
      ini({
        tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .fit <- .nlmixr(one.cmt, theo_sd, "agq", agqControl(print = 0, maxOuterIterations = 0L, covMethod = ""))
    expect_equal(.fit$control$nAGQ, 2L)
    .ll <- .fit$llikObs
    .ll <- .ll[!is.na(.ll)] # dose records
    # llikObs of a Gaussian row omits the 2*pi constant of dnorm()
    expect_equal(
      .ll,
      dnorm(.fit$DV, .fit$IPRED, .fit$theta[["add.sd"]], log = TRUE) + 0.5 * log(2 * pi),
      tolerance = 1e-10
    )
  })

  test_that("forward-difference S scores are taken from the estimates", {
    # covDerivMethod = "forward" differences each subject's -2LL from its value at
    # the estimates.  That value was read after the pooled gradient's own
    # forward legs, so it was the last leg's, and every score carried that leg's
    # shift: S came out 100 to 1000 times the central-difference S.
    one.cmt <- function() {
      ini({
        tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .ctl <- function(...) {
      foceiControl(print = 0, maxOuterIterations = 0L, covMethod = "s", covFull = FALSE, ...)
    }
    .central <- .nlmixr(one.cmt, theo_sd, "focei", .ctl())
    .forward <- .nlmixr(one.cmt, theo_sd, "focei", .ctl(covDerivMethod = "forward"))
    expect_equal(.forward$covMethod, "s")
    # forward differences agree with central ones to their truncation error
    expect_equal(diag(.forward$env$S0), diag(.central$env$S0), tolerance = 0.25)
  })

  test_that("covDerivMethod = \"forward\" keeps the R matrix", {
    # the R matrix is always the central stencil; foceiCalcR() stopped on a
    # forward derivative method instead, so "r" and "r,s" failed
    one.cmt <- function() {
      ini({
        tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .ctl <- function(...) {
      foceiControl(print = 0, maxOuterIterations = 0L, covFull = FALSE, ...)
    }
    .central <- .nlmixr(one.cmt, theo_sd, "focei", .ctl())
    .forward <- .nlmixr(one.cmt, theo_sd, "focei", .ctl(covDerivMethod = "forward"))
    expect_equal(.central$covMethod, "r,s")
    expect_equal(.forward$covMethod, "r,s")
    expect_identical(.forward$env$R.0, .central$env$R.0)
    # the sandwich differs only through the forward-difference S
    expect_equal(sqrt(diag(.forward$cov)), sqrt(diag(.central$cov)), tolerance = 0.2)
    .r <- .nlmixr(one.cmt, theo_sd, "focei", .ctl(covMethod = "r", covDerivMethod = "forward"))
    expect_equal(.r$covMethod, "r")
    expect_equal(unname(.r$cov), unname(.central$covR))
  })

  test_that("a requested full sandwich gets its S when the native step falls back to R", {
    # One subject: the theta-only S is rank one, so the native "r,s" falls back to
    # the R matrix.  The full covariance then computed no S, and the requested
    # "r,s (full)" was dropped without a word; it is now computed, checked, and
    # the reason it is not installed is reported.
    one.cmt <- function() {
      ini({
        tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .f <- .nlmixr(
      one.cmt,
      theo_sd[theo_sd$ID == 1, ],
      "focei",
      foceiControl(print = 0, maxOuterIterations = 0L, cholAccept = 0)
    )
    expect_equal(.covFdType(.f$covMethod), "r")
    expect_false(.covIsFull(.f$covMethod))
    expect_true(is.matrix(.f$env$.fdFullS))
    expect_true(any(startsWith(.f$runInfo, "\"r,s (full)\" covariance ")))
  })

  test_that("an all-zero R is not corrected into a covariance", {
    # the objective does not move with tz at all (its covariate is 0), so R = 0;
    # cholSE0 adds cholSEtol to a zero 1x1 with nothing to scale it by, which
    # passed as a small "r+" correction and installed 1/cholSEtol
    d <- data.frame(ID = rep(1:2, each = 3), TIME = rep(1:3, 2), DV = 5, Z = 0)
    flat <- function() {
      ini({
        ta <- fix(5)
        tz <- 0.5
        add.sd <- fix(1)
      })
      model({
        cp <- ta + tz * Z
        cp ~ add(add.sd)
      })
    }
    .f <- .nlmixr(flat, d, "focei", foceiControl(print = 0, maxOuterIterations = 0L))
    expect_equal(.f$env$R.0[1, 1], 0)
    expect_equal(.f$covMethod, "failed")
    expect_null(.f$cov)
  })

  test_that("an all-zero R of several thetas is not corrected into a covariance", {
    # the zero diagonal leaves cholSE0 nothing to scale its correction by for any
    # theta, not only a single one
    d <- data.frame(ID = rep(1:2, each = 3), TIME = rep(1:3, 2), DV = 5, Z = 0, W = 0)
    flat2 <- function() {
      ini({
        ta <- fix(5)
        tz <- 0.5
        tw <- 0.3
        add.sd <- fix(1)
      })
      model({
        cp <- ta + tz * Z + tw * W
        cp ~ add(add.sd)
      })
    }
    .f <- .nlmixr(flat2, d, "focei", foceiControl(print = 0, maxOuterIterations = 0L))
    expect_equal(unname(.f$env$R.0), matrix(0, 2, 2))
    expect_equal(.f$covMethod, "failed")
    expect_null(.f$cov)
  })

  test_that("a non-positive-definite R or S is never installed as it is", {
    # one estimated parameter at a point where the objective is concave: R < 0,
    # which cholSE0 (like for every 1x1 matrix) called positive definite, so
    # 1/(cholSEtol*|R|) was installed as "r"
    d <- data.frame(ID = rep(1:2, each = 3), TIME = rep(1:3, 2), DV = 5)
    peak <- function() {
      ini({
        ta <- 3
        add.sd <- fix(1)
      })
      model({
        cp <- 10 * exp(-(ta - 3)^2)
        cp ~ add(add.sd)
      })
    }
    .f1 <- .nlmixr(peak, d, "focei", foceiControl(print = 0, maxOuterIterations = 0L))
    expect_lt(.f1$env$R.0[1, 1], 0)
    expect_equal(.f1$covMethod, "|r|")
    expect_equal(.f1$cov[1, 1], 1 / abs(.f1$env$R.0[1, 1]))
    # one subject: S is rank one and cannot be repaired; it was labelled "s"
    # with no covariance
    one.cmt <- function() {
      ini({
        tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .f2 <- .nlmixr(
      one.cmt,
      theo_sd[theo_sd$ID == 1, ],
      "focei",
      foceiControl(print = 0, maxOuterIterations = 0L, covMethod = "s", cholAccept = 0)
    )
    expect_true(!is.null(.f2$cov) || identical(.f2$covMethod, "failed"))
    expect_equal(.f2$covMethod, "failed")
  })

  # theo_sd one-compartment ODE model, with the covariance solves (covSolveTol) and the
  # inner problem (trustFterm/trustMterm) tight: the marginal objective is then smooth to
  # ~1e-6 along a finite-difference step.  At the default tolerances (ODE rtol = 1e-3) its
  # noise is ~1e-3, as large as what a step changes it by, and no finite-difference
  # Hessian can be checked against anything.
  .quietOneCmt <- function() {
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
  .quietCtl <- function(...) {
    foceiControl(print = 0, covSolveTol = 1e-9, trustFterm = 1e-8, trustMterm = 1e-8, ...)
  }

  test_that("the finite-difference R matrix is the Hessian of the marginal objective", {
    skip_on_cran()
    # Every leg optimizes the ETAs again, starting from the ones the fit converged to, and
    # the centre of the stencil is the objective evaluated the same way, so R must be the
    # analytic observed information.  The same finite differences move by up to 4% when
    # the step is made 4 times larger or the tolerances 100 times tighter; 5% is the bound.
    # A centre taken from the final objective (ETAs optimized by another procedure, at the
    # fit's ODE tolerance) makes every diagonal of R negative here, at the initial
    # estimates
    .an0 <- .nlmixr(.quietOneCmt, theo_sd, "focei", .quietCtl(maxOuterIterations = 0L, covMethod = "analytic"))
    .r0 <- .nlmixr(.quietOneCmt, theo_sd, "focei", .quietCtl(maxOuterIterations = 0L, covMethod = "r", covFull = FALSE))
    .th <- c("tka", "tcl", "tv", "add.sd")
    # Omega is held at its value in R: the theta block of the information
    .info <- solve(.an0$cov)[.th, .th]
    expect_identical(.r0$covMethod, "r")
    expect_lt(max(abs(diag(.r0$env$R.0) / diag(.info) - 1)), 0.05)
    # the full (covFull) stage is the same derivative over theta, sigma and Omega; at the
    # estimates, where it is positive definite
    .an <- .nlmixr(.quietOneCmt, theo_sd, "focei", .quietCtl(covMethod = "analytic"))
    .rf <- .nlmixr(.quietOneCmt, theo_sd, "focei", .quietCtl(covMethod = "r"))
    expect_identical(.rf$covMethod, "r (full)")
    .nm <- rownames(.an$cov)
    .rel <- sqrt(diag(.rf$cov))[.nm] / sqrt(diag(.an$cov)) - 1
    # om.eta.ka has no step search in the full shape; its SE is 1% off here and 12% on
    # some CI platforms (issue #1140, A6), so it gets its own bound
    for (.p in .nm) {
      .bound <- if (.p == "om.eta.ka") 0.15 else 0.05
      expect_lt(abs(.rel[[.p]]), .bound, label = sprintf("%s SE relative error", .p))
    }
  })

  test_that("a finite-difference covariance of ETAs held fixed says it is conditional on them", {
    one.cmt <- function() {
      ini({
        tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .eta <- .nlmixr(one.cmt, theo_sd, "focei", foceiControl(print = 0, covMethod = ""))$etaMat
    .msg <- "covariance is conditional on the ETAs; setCov() is marginal"
    # maxInnerIterations = 0 evaluates the given ETAs, in the fit and in its covariance
    .fixed <- .nlmixr(
      one.cmt,
      theo_sd,
      "posthoc",
      posthocControl(
        print = 0,
        maxInnerIterations = 0L,
        etaMat = .eta,
        covMethod = "r",
        covFull = FALSE
      )
    )
    expect_identical(.covFdType(.fixed$covMethod), "r")
    expect_true(any(grepl(.msg, .fixed$runInfo, fixed = TRUE)))
    # optimizing them: the covariance of the marginal likelihood, nothing to say
    .opt <- .nlmixr(
      one.cmt,
      theo_sd,
      "posthoc",
      posthocControl(
        print = 0,
        maxInnerIterations = 1000L,
        etaMat = .eta,
        covMethod = "r",
        covFull = FALSE
      )
    )
    expect_identical(.covFdType(.opt$covMethod), "r")
    expect_false(any(grepl(.msg, .opt$runInfo, fixed = TRUE)))
    # recomputing it (getVarCov(force = TRUE), as setCov()) gives the marginal one
    .v <- suppressMessages(suppressWarnings(nlme::getVarCov(.fixed, force = TRUE)))
    expect_equal(sqrt(diag(.v)), sqrt(diag(.opt$cov)), tolerance = 1e-6)
  })

  test_that("at the default tolerances the finite-difference R is the analytic information", {
    skip_on_cran()
    # The probes run at the fit's tolerances times 1e-3 (ODE rtol 1e-7, atol 1e-9, inner
    # 1e-9 at sigdig 3), where every SE is within 4% of the analytic one; at the fit's own
    # (rtol 1e-3, inner 1e-5) the probes difference numerical noise
    .an <- .nlmixr(.quietOneCmt, theo_sd, "focei", foceiControl(print = 0, covMethod = "analytic"))
    .r <- .nlmixr(.quietOneCmt, theo_sd, "focei", foceiControl(print = 0, covMethod = "r", covFull = FALSE))
    expect_identical(.r$covMethod, "r")
    .th <- c("tka", "tcl", "tv", "add.sd")
    # Omega is held at its value in a theta-only R: the theta block of the information
    .seCond <- sqrt(diag(solve(solve(.an$cov)[.th, .th])))
    expect_lt(max(abs(sqrt(diag(.r$cov))[.th] / .seCond - 1)), 0.05)
  })

  test_that("the covariance step tightens the inner problem of every inner optimizer", {
    skip_on_cran()
    # each inner optimizer reads its own convergence tolerances (n1qn1 epsilon,
    # lbfgsb3c's pgtol/abstol/reltol/factr); the step tightens all of them, so the
    # finite-difference R is the analytic information whichever one runs the legs
    # (measured 0.7% for n1qn1 and 2.8% for lbfgsb3c, against 32% and 66% with the
    # inner problems at the estimation's tolerances)
    .th <- c("tka", "tcl", "tv", "add.sd")
    for (.io in c("n1qn1", "lbfgsb3c")) {
      .an <- .nlmixr(.quietOneCmt, theo_sd, "focei", foceiControl(print = 0, innerOpt = .io, covMethod = "analytic"))
      .r <- .nlmixr(
        .quietOneCmt,
        theo_sd,
        "focei",
        foceiControl(print = 0, innerOpt = .io, covMethod = "r", covFull = FALSE)
      )
      expect_identical(.r$covMethod, "r", label = .io)
      .rel <- sqrt(diag(.r$cov))[.th] / sqrt(diag(.an$cov))[.th] - 1
      expect_lt(max(abs(.rel)), 0.05, label = .io)
    }
  })

  test_that("the covariance probe tolerances of a high-sigdig fit stay solvable", {
    skip_on_cran()
    # sigdig = 10 asks for ODE and inner tolerances near 1e-11; 1e-3 of those is floored
    # at 1e-14
    .f <- .nlmixr(
      .quietOneCmt,
      theo_sd,
      "focei",
      foceiControl(print = 0, sigdig = 10, covMethod = "r", covFull = FALSE, maxOuterIterations = 2L)
    )
    expect_identical(.f$covMethod, "r")
    expect_true(all(is.finite(sqrt(diag(.f$cov)))))
  })

  test_that("the covariance step runs at its probe tolerances and leaves estimation as it was", {
    skip_on_cran()
    # Estimation, its objective, ETAs and tables are those of the fit's own tolerances:
    # the covariance step that follows tightens its solves and inner problems only for
    # itself.
    .none <- .nlmixr(.quietOneCmt, theo_sd, "focei", foceiControl(print = 0, covMethod = ""))
    .cov <- .nlmixr(.quietOneCmt, theo_sd, "focei", foceiControl(print = 0))
    expect_identical(.cov$objf, .none$objf)
    expect_identical(.cov$theta, .none$theta)
    expect_identical(.cov$omega, .none$omega)
    expect_identical(.cov$etaMat, .none$etaMat)
    expect_identical(as.data.frame(.cov)$IPRED, as.data.frame(.none)$IPRED)
    expect_identical(as.data.frame(.cov)$CWRES, as.data.frame(.none)$CWRES)
    # NULL covSolveTol is the fit's atol and rtol times 1e-3, capped at 1e-7: a fit at
    # atol = rtol = 1e-4 gets 1e-7 for both, as covSolveTol = 1e-7 sets them
    .ctl <- function(...) {
      foceiControl(
        print = 0,
        maxOuterIterations = 0L,
        covMethod = "r",
        covFull = FALSE,
        rxControl = rxode2::rxControl(atol = 1e-4, rtol = 1e-4),
        ...
      )
    }
    .rule <- .nlmixr(.quietOneCmt, theo_sd, "focei", .ctl())
    expect_identical(.rule$cov, .nlmixr(.quietOneCmt, theo_sd, "focei", .ctl(covSolveTol = 1e-7))$cov)
    # the inner tolerance is the fit's trustFterm/trustMterm (1e-5 at sigdig 3) times
    # 1e-3, capped at 1e-9
    expect_identical(.rule$cov, .nlmixr(.quietOneCmt, theo_sd, "focei", .ctl(covInnerTol = 1e-9))$cov)
    # and each setting is used: another value gives another matrix
    expect_false(identical(.rule$cov, .nlmixr(.quietOneCmt, theo_sd, "focei", .ctl(covSolveTol = 1e-6))$cov))
    expect_false(identical(.rule$cov, .nlmixr(.quietOneCmt, theo_sd, "focei", .ctl(covInnerTol = 1e-7))$cov))
  })
})
