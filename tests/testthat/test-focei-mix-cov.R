nmTest({

  ## Standard errors for mixture proportions -- NONMEM 7 Technical Guide
  ## eq. (7.51)-(7.54), which build the mixture parameters' information matrix
  ## from the same per-subject scores g_ia as the gradient (1.194)/(1.197).
  ##
  ## Two defects are guarded here.  The proportions were forced out of the
  ## covariance entirely (skipCov), so they reported SE = NA under every method.
  ## And foceiS() built each subject's score by differencing component 0's
  ## likelihood rather than the marginal, which made the S matrix singular for
  ## EVERY mixture model -- a parameter that only enters another component got an
  ## exactly-zero score, and covMethod="r,s" silently degraded to "r".

  test_that("the mixture covariance block is rotated by the FULL mexpit Jacobian", {
    ## No fit: the rotation is pure linear algebra, and the full Jacobian (not
    ## its diagonal) is what carries the proportions' cross-covariances onto the
    ## reported scale.
    .t <- c(0.3, -0.8)
    .p <- rxode2::mexpit(.t)
    ## numeric Jacobian of mexpit at .t, for comparison with diag(p) - p p'
    .num <- vapply(seq_along(.t), function(j) {
      .h <- 1e-6
      .tp <- .t; .tp[j] <- .t[j] + .h
      .tm <- .t; .tm[j] <- .t[j] - .h
      (rxode2::mexpit(.tp) - rxode2::mexpit(.tm)) / (2 * .h)
    }, numeric(length(.t)))
    .J <- diag(.p, nrow = length(.p)) - outer(.p, .p)
    expect_equal(.J, .num, tolerance = 1e-6)

    ## and the rotation is A cov A' with A the identity off the mixture block
    .nm <- c("tcl", "p1", "p2")
    .cov <- matrix(c(4, 1, 2,
                     1, 9, 3,
                     2, 3, 16), 3, 3, dimnames = list(.nm, .nm))
    .out <- .mixCovToProbScale(.cov, c("p1", "p2"), .p)
    .A <- diag(1, 3)
    .A[2:3, 2:3] <- .J
    expect_equal(unname(.out), unname(.A %*% .cov %*% t(.A)))
    expect_equal(dimnames(.out), dimnames(.cov))
    ## a non-mixture cov, or names that are not present, is left alone
    expect_equal(.mixCovToProbScale(.cov, character(0), numeric(0)), .cov)
    expect_equal(.mixCovToProbScale(.cov, "nope", .p[1]), .cov)
  })

  ## Two well-separated components with everything but the proportion and the
  ## residual sd fixed at truth.  The responsibilities are then essentially 0/1,
  ## so p1 is the observed group fraction and its variance is p(1-p)/n -- an
  ## analytic target the reported SE can be held to without any replicate study.
  .mixCovData <- function(nSub = 60L, pTrue = 0.45, clTrue = c(1.0, 8.0), seed = 1001L) {
    set.seed(seed)
    .grp <- sample.int(2L, nSub, TRUE, prob = c(pTrue, 1 - pTrue))
    .sim <- rxode2::rxode2({
      ka <- 1.1
      cl <- CLI
      v <- 20
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
    })
    .ev <- rxode2::et(rxode2::et(amt = 320, cmt = "depot"),
                      c(0.25, 0.5, 1, 2, 4, 8, 12, 24))
    .obs <- do.call(rbind, lapply(seq_len(nSub), function(i) {
      .s <- rxode2::rxSolve(.sim, params = c(CLI = clTrue[.grp[i]]), .ev,
                            returnType = "data.frame")
      data.frame(ID = i, TIME = .s$time,
                 DV = .s$cp + stats::rnorm(nrow(.s), 0, 0.05), AMT = 0, EVID = 0)
    }))
    .d <- rbind(data.frame(ID = seq_len(nSub), TIME = 0, DV = NA_real_,
                           AMT = 320, EVID = 1), .obs)
    list(data = .d[order(.d$ID, .d$TIME, -.d$EVID), ], n = nSub)
  }

  .mixCovMod <- function() {
    ini({
      tka <- fix(log(1.1))
      tcl1 <- fix(log(1.0))
      tcl2 <- fix(log(8.0))
      tv <- fix(log(20))
      p1 <- 0.45
      eta.cl ~ fix(0.01)
      add.sd <- 0.05
    })
    model({
      ka <- exp(tka)
      cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("every focei covariance method gives the mixture proportion a calibrated SE", {
    .dat <- .mixCovData()
    .fits <- lapply(c("r,s", "r", "s"), function(.cm) {
      suppressWarnings(nlmixr2(.mixCovMod, .dat$data, "focei",
        foceiControl(print = 0, outerOpt = "lbfgsb3c", maxOuterIterations = 200L,
                     maxInnerIterations = 100L, covMethod = .cm, calcTables = FALSE)))
    })
    for (.f in .fits) {
      .p <- .f$parFixedDf["p1", "Estimate"]
      .se <- .f$parFixedDf["p1", "SE"]
      ## reported at all -- these were NA for every mixture fit
      expect_true(is.finite(.se) && .se > 0)
      ## and on the PROBABILITY scale the estimate is reported on: the mlogit
      ## SE here is ~0.26, four times larger, so this pins the Jacobian rotation
      expect_equal(unname(.se), sqrt(.p * (1 - .p) / .dat$n), tolerance = 0.05)
      ## the proportion is in $cov by name
      expect_true("p1" %in% rownames(.f$cov))
    }
    ## the three methods agree with each other
    .ses <- vapply(.fits, function(.f) .f$parFixedDf["p1", "SE"], numeric(1))
    expect_equal(max(.ses) / min(.ses), 1, tolerance = 0.05)
  })

  test_that("the rotation is keyed on the theta slot, not on ui$mixProbs order", {
    ## ui$mixProbs follows the mix() CALL; thetaMixIndex follows ini().  When the
    ## two orders differ, everything downstream -- the covariance's rows,
    ## op_focei.mixProb, $mixProbabilities, se/popDf -- is keyed on the THETA
    ## slot, so indexing a theta-ordered covariance by mixProbs puts the Jacobian
    ## on the wrong rows.
    .mod <- function() {
      ini({
        tka <- log(1.1)
        p2 <- 0.30          # deliberately declared BEFORE p1
        p1 <- 0.60
        tcl1 <- log(1)
        tcl2 <- log(8)
        tcl3 <- log(30)
        tv <- log(20)
        eta.cl ~ 0.01
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl), p2,
                  exp(tcl3 + eta.cl))
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::assertRxUi(.mod))
    .slot <- names(.ui$theta)[.ui$thetaMixIndex]
    ## thetaMixIndex is in COMPONENT order, so the slots it names are mixProbs
    expect_equal(.slot, .ui$mixProbs)
    ## ...but this model's ini() declares them in the OTHER order, so those slot
    ## POSITIONS are descending.  That is what makes the test non-trivial: the
    ## covariance's rows are built in ascending theta order (p2 then p1) while
    ## the Jacobian is built in component order (p1 then p2), so the two must be
    ## bridged by name rather than by position.
    expect_gt(.ui$thetaMixIndex[1], .ui$thetaMixIndex[2])

    ## component-ordered probabilities, far enough apart that a swap shows
    .p <- c(0.60, 0.30)
    .e <- new.env(parent = emptyenv())
    .e$ui <- .ui
    .e$mixProbabilities <- c(.p, 1 - sum(.p))
    ## the covariance is in THETA order, so its mixture rows are named p2, p1
    .nm <- c(.slot, "add.sd")
    .e$cov <- diag(c(1, 4, 9))
    dimnames(.e$cov) <- list(.nm, .nm)
    .mixInstallProbScaleCov(.e)

    ## component m's Jacobian must land on the row for theta slot mixIdx[m]
    .j <- diag(.p, nrow = 2L) - outer(.p, .p)
    .a <- diag(1, 3)
    .a[1:2, 1:2] <- .j
    .want <- .a %*% diag(c(1, 4, 9)) %*% t(.a)
    expect_equal(unname(.e$cov), unname(.want))
    ## and concretely: the two mixture variances are NOT interchangeable here
    expect_false(isTRUE(all.equal(.e$cov[1, 1], .e$cov[2, 2])))
  })

  test_that("a fix()ed proportion is not given an appended variance", {
    ## thetaMixIndex still lists a fix()ed proportion, but it was never
    ## estimated -- appending a row for it would report a non-zero SE for a
    ## parameter with no uncertainty.
    .m <- function() {
      ini({
        tcl <- log(1)
        p1 <- fix(0.40)
        p2 <- 0.35
        tv <- log(20)
        eta.cl ~ 0.01
        add.sd <- 0.05
      })
      model({
        cl <- mix(exp(tcl + eta.cl), p1, exp(tcl + 2 + eta.cl), p2,
                  exp(tcl + 3 + eta.cl))
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::assertRxUi(.m))
    .n <- 100L
    .k <- round(.n * 0.35)
    .r <- cbind(rep(0.40, .n),
                c(rep(1, .k), rep(0, .n - .k)))
    .r <- cbind(.r, 1 - rowSums(.r))
    .e <- new.env(parent = emptyenv())
    .e$ui <- .ui
    .e$mixProbabilities <- c(0.40, .k / .n, 1 - 0.40 - .k / .n)
    .nm <- c("tcl", "add.sd")
    .e$cov <- diag(c(4, 9)); dimnames(.e$cov) <- list(.nm, .nm)
    .e$mixList <- lapply(seq_len(3L), function(.j) data.frame(prob = .r[, .j]))
    .mixCovAppendBlock(.e)
    ## exactly one appended row, for the ESTIMATED proportion
    expect_equal(rownames(.e$cov), c("tcl", "add.sd", "p2"))
    expect_false("p1" %in% rownames(.e$cov))
    expect_equal(unname(sqrt(diag(.e$cov))[3]),
                 sqrt(0.35 * 0.65 / .n), tolerance = 1e-6)
  })

  test_that("a mixture proportion's CI stays inside (0, 1)", {
    ## The estimate IS a probability, so the generic symmetric
    ## backTransform(est +/- z*SE) interval walks straight out of (0, 1) -- a
    ## real fit reported p1 = 0.648 (-0.045, 1.34).  It was also built from the
    ## covariance BEFORE the probability-scale rotation, so it did not even
    ## agree with the SE printed beside it.
    .ui <- rxode2::rxUiDecompress(rxode2::assertRxUi(.mixCovMod))
    .pf <- data.frame(Estimate = c(0.648, 0.5), SE = c(0.0862, 0.30),
                      `CI Lower` = c(-0.045, -0.1), `CI Upper` = c(1.34, 1.1),
                      row.names = c("p1", "add.sd"), check.names = FALSE)
    .out <- .mixParFixedCi(.ui, .pf, 0.95)
    expect_true(.out["p1", "CI Lower"] > 0 && .out["p1", "CI Upper"] < 1)
    ## it is the logit-scale interval, built from the REPORTED SE
    .p <- 0.648; .s <- 0.0862; .j <- .p * (1 - .p)
    expect_equal(.out["p1", "CI Lower"],
                 rxode2::expit(rxode2::logit(.p) - 1.959964 * .s / .j),
                 tolerance = 1e-5)
    ## asymmetric about the estimate, which is the honest shape here
    expect_false(isTRUE(all.equal(.p - .out["p1", "CI Lower"],
                                  .out["p1", "CI Upper"] - .p)))
    ## a non-mixture row is untouched
    expect_equal(.out["add.sd", "CI Lower"], -0.1)
  })

  test_that("a proportion at the boundary is flagged, not reported as precise", {
    ## The probability-scale SE carries a factor p(1-p), so it goes to ZERO as a
    ## proportion approaches 0 or 1.  That is the correct delta-method answer but
    ## reads as certainty, so it must be said out loud.
    expect_warning(.mixWarnBoundary(c(0.001, 0.999)), "near 0/1")
    expect_warning(.mixWarnBoundary(c(0.45, 0.996)), "near 0/1")
    expect_silent(.mixWarnBoundary(c(0.45, 0.55)))
    expect_silent(.mixWarnBoundary(c(0.02, 0.98)))
    expect_silent(.mixWarnBoundary(numeric(0)))
    expect_silent(.mixWarnBoundary(c(NA_real_, 0.5)))
  })

  test_that("a proportion missing from the covariance does not block the others", {
    ## A fix()ed proportion is dropped from the covariance by skipCov.  Bailing
    ## on the whole rotation when a name is absent left the REMAINING
    ## proportions reported on the mlogit scale -- measured 0.265 where the
    ## probability scale is 0.063, a factor of 1/(p(1-p)).
    .nm <- c("p2", "add.sd")
    .cov <- diag(c(4, 9)); dimnames(.cov) <- list(.nm, .nm)
    .p <- c(0.40, 0.35)                       # p1 (fixed, absent) and p2
    .out <- .mixCovToProbScale(.cov, c("p1", "p2"), .p)
    ## only p2's row is rotated, by its own diagonal Jacobian entry p2(1-p2)
    .j22 <- .p[2] * (1 - .p[2])
    expect_equal(unname(.out[1, 1]), .j22^2 * 4)
    expect_equal(unname(.out[2, 2]), 9)       # add.sd untouched
    ## and it is NOT left unrotated
    expect_false(isTRUE(all.equal(unname(.out[1, 1]), 4)))

    ## with every name absent the matrix is returned as-is
    .nm2 <- c("tcl", "add.sd")
    .cov2 <- diag(c(4, 9)); dimnames(.cov2) <- list(.nm2, .nm2)
    expect_equal(.mixCovToProbScale(.cov2, c("p1", "p2"), .p), .cov2)
  })

  test_that("setCov() round trips without re-rotating the mixture block", {
    ## setCov() re-installs a CACHED covariance (covList) by handing it back as a
    ## matrix, which refits and would rotate an already-probability-scale matrix
    ## a second time -- shrinking the proportion's SE by p(1-p) every round trip
    ## (measured 0.0644 -> 0.0155).
    .dat <- .mixCovData()
    .f <- suppressWarnings(nlmixr2(.mixCovMod, .dat$data, "focei",
      foceiControl(print = 0, outerOpt = "lbfgsb3c", maxOuterIterations = 200L,
                   maxInnerIterations = 100L, covMethod = "r,s", calcTables = FALSE)))
    .se0 <- .f$parFixedDf["p1", "SE"]
    ## without this the round-trip check passes vacuously against the old code,
    ## where the SE was NA at both ends and NA == NA
    expect_true(is.finite(.se0) && .se0 > 0)
    setCov(.f, "s")
    expect_true("r,s" %in% names(.f$env$covList))
    setCov(.f, "r,s")                       # served from covList, not recomputed
    .se1 <- .f$parFixedDf["p1", "SE"]
    expect_true(is.finite(.se1) && .se1 > 0)
    expect_equal(unname(.se1), unname(.se0), tolerance = 1e-10)
  })

  test_that("the S matrix is no longer singular for a mixture model", {
    .dat <- .mixCovData()
    .f <- suppressWarnings(nlmixr2(.mixCovMod, .dat$data, "focei",
      foceiControl(print = 0, outerOpt = "lbfgsb3c", maxOuterIterations = 200L,
                   maxInnerIterations = 100L, covMethod = "r,s", calcTables = FALSE)))
    ## "r,s" is honoured rather than degraded to "r" by a non-PD S
    expect_true(isTRUE(.f$env$S.pd))
    expect_equal(.f$covMethod, "r,s")
    ## no structurally-zero score direction
    expect_true(all(diag(.f$env$S0) > 1e-6))

    ## and the mixture block IS the outer product of the analytic per-subject
    ## score -- NONMEM (7.51).  S = 0.25 * sum_i g_i g_i' with
    ## g_i[l] = -2*(r_il - pi_l), so the block is sum_i (r_il - pi_l)(r_im - pi_m).
    ##
    ## Relative tolerance, not exact: $mixList's responsibilities are the FINAL
    ## table's, while S's were taken at the covariance step's base point, and the
    ## two states differ slightly on a converged fit.  1e-2 still pins the
    ## identity -- differencing component 0 alone (the defect) puts this entry
    ## orders of magnitude out, not fractions of a percent.
    .pi <- .f$env$mixProbabilities
    .R <- do.call(cbind, lapply(.f$env$mixList, function(z) z$prob))
    .free <- seq_len(length(.pi) - 1L)
    .D <- sweep(.R[, .free, drop = FALSE], 2, .pi[.free], "-")
    .i <- match(.f$ui$mixProbs, rownames(.f$cov))
    expect_equal(unname(.f$env$S0[.i, .i, drop = FALSE]), unname(t(.D) %*% .D),
                 tolerance = 1e-2)
  })

  test_that("covMethod='imp' gives the mixture proportion the same calibrated SE", {
    .dat <- .mixCovData()
    .f <- suppressWarnings(nlmixr2(.mixCovMod, .dat$data, "imp",
      impmapControl(print = 0, nIter = 8L, covMethod = "imp", calcTables = FALSE)))
    expect_equal(.f$covMethod, "imp")
    .p <- .f$parFixedDf["p1", "Estimate"]
    .se <- .f$parFixedDf["p1", "SE"]
    ## the imp MC covariance built its proposals from component 0 only, so the
    ## proportion's direction was exactly flat and this came back 0
    expect_true(is.finite(.se) && .se > 0)
    expect_equal(unname(.se), sqrt(.p * (1 - .p) / .dat$n), tolerance = 0.05)
  })

  ## The appended block for engines whose covariance has no mixture rows at all
  ## (saem excludes the proportions from its kernel parameter vector).  Driven
  ## from a synthetic env so both branches are exercised without a fit.
  .mkMixEnv <- function(r1, pi1, n = 100L) {
    ## a REAL ui: the covariance helpers key the proportions on their theta slot
    ## (names(ui$theta)[ui$thetaMixIndex]), so a stub list will not do
    .m <- function() {
      ini({
        tcl <- log(1)
        p1 <- 0.45
        tv <- log(20)
        eta.cl ~ 0.01
        add.sd <- 0.05
      })
      model({
        cl <- mix(exp(tcl + eta.cl), p1, exp(tcl + 2 + eta.cl))
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .e <- new.env(parent = emptyenv())
    .e$ui <- rxode2::rxUiDecompress(rxode2::assertRxUi(.m))
    .nm <- c("tcl", "tv")
    .e$cov <- matrix(c(4, 1, 1, 9), 2, 2, dimnames = list(.nm, .nm))
    .e$mixProbabilities <- c(pi1, 1 - pi1)
    .e$mixList <- list(data.frame(prob = r1), data.frame(prob = 1 - r1))
    .e
  }

  test_that("a covariance with no mixture rows gets the NONMEM (7.51) block appended", {
    ## hard 0/1 responsibilities with mean == the proportion: the EM fixed point,
    ## where the information is n*p*(1-p) and the answer is the binomial variance
    .n <- 100L
    .r <- c(rep(1, 45), rep(0, 55))
    .e <- .mkMixEnv(.r, 0.45, .n)
    .mixCovAppendBlock(.e)
    expect_true("p1" %in% rownames(.e$cov))
    expect_equal(dim(.e$cov), c(3L, 3L))
    expect_equal(unname(sqrt(diag(.e$cov))[3]), sqrt(0.45 * 0.55 / .n), tolerance = 1e-8)
    ## the pre-existing block is untouched and the new one is uncorrelated with
    ## it -- the cross terms (7.52)-(7.54) need per-subject scores this engine
    ## does not expose, so they are deliberately absent
    expect_equal(unname(.e$cov[1:2, 1:2]), matrix(c(4, 1, 1, 9), 2, 2))
    expect_equal(unname(.e$cov[1:2, 3]), c(0, 0))
  })

  test_that("the appended block is refused when the fit is not at the EM fixed point", {
    ## p != mean_i r_i: an information matrix reports the precision of an MLE,
    ## and this is not one.  saem lands here (nlmixr2est#1058), and reporting a
    ## number would be a confident-looking SE on a non-stationary estimate.
    .e <- .mkMixEnv(c(rep(1, 30), rep(0, 70)), 0.64, 100L)
    expect_warning(.mixCovAppendBlock(.e), "score-zero")
    expect_false("p1" %in% rownames(.e$cov))
    expect_equal(dim(.e$cov), c(2L, 2L))
  })

  test_that("the stationarity gate does not loosen with the number of subjects", {
    ## The gate is the score statistic s' I^-1 s, not |mean(r) - p|: the score is
    ## N*(mean(r) - p), so an absolute tolerance on the mean accepts a score of 1
    ## at N=100 and 100 at N=10000.  Hold the DEVIATION fixed and grow N; the
    ## same deviation must be refused at least as firmly at the larger N.
    .dev <- 0.02
    .mk <- function(n) {
      .k <- round(n * (0.45 + .dev))
      .mkMixEnv(c(rep(1, .k), rep(0, n - .k)), 0.45, n)
    }
    for (.n in c(100L, 2000L)) {
      .e <- .mk(.n)
      expect_warning(.mixCovAppendBlock(.e), "score-zero")
      expect_false("p1" %in% rownames(.e$cov))
    }
    ## and a fit that IS at the fixed point is still accepted at both sizes
    for (.n in c(100L, 2000L)) {
      .k <- round(.n * 0.45)
      .e <- .mkMixEnv(c(rep(1, .k), rep(0, .n - .k)), .k / .n, .n)
      .mixCovAppendBlock(.e)
      expect_true("p1" %in% rownames(.e$cov))
      expect_equal(unname(sqrt(diag(.e$cov))[3]),
                   sqrt(0.45 * 0.55 / .n), tolerance = 1e-6)
    }
  })

  test_that("covMethod='analytic' declines a mixture rather than reporting one component", {
    ## The augmented sensitivity model differentiates ONE component's conditional
    ## likelihood, not the marginal, and has no mixture-proportion block at all.
    .dat <- .mixCovData()
    .f <- suppressWarnings(nlmixr2(.mixCovMod, .dat$data, "focei",
      foceiControl(print = 0, outerOpt = "lbfgsb3c", maxOuterIterations = 200L,
                   maxInnerIterations = 100L, covMethod = "analytic",
                   calcTables = FALSE)))
    expect_false(identical(.covBaseName(.f$covMethod), "analytic"))
    ## and the fallback still reports the proportion
    expect_true(is.finite(.f$parFixedDf["p1", "SE"]))
  })

})
