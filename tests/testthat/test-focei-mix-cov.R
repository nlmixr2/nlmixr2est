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
    expect_equal(.mixCovToProbScale(.cov, c("nope"), .p[1]), .cov)
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

  test_that("covMethod='analytic' declines a mixture rather than reporting one component", {
    ## The augmented sensitivity model differentiates ONE component's conditional
    ## likelihood, not the marginal, and has no mixture-proportion block at all.
    .dat <- .mixCovData()
    .f <- suppressWarnings(nlmixr2(.mixCovMod, .dat$data, "focei",
      foceiControl(print = 0, outerOpt = "lbfgsb3c", maxOuterIterations = 200L,
                   maxInnerIterations = 100L, covMethod = "analytic",
                   calcTables = FALSE)))
    expect_false(identical(.f$covMethod, "analytic"))
    ## and the fallback still reports the proportion
    expect_true(is.finite(.f$parFixedDf["p1", "SE"]))
  })

})
