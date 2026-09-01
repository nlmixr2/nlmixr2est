# Fit-path regression test for rxode2's linCmt() transition-matrix engage
# rule.  The rule assembles a matrix only on evidence that a row's interval
# RECURS in the design, so a design whose gaps never repeat must build none.
#
# Why this test lives here rather than in rxode2: the rule once counted any
# delta-memo hit as evidence of recurrence, and in a fit a single row reaches
# the row-Jacobian several times (the generated value line runs from dydt and
# from calc_lhs, and the inner problem re-walks each subject), so a row asking
# twice looked like a recurrence -- a strictly non-repeating design built
# millions of matrices in a fit while building zero in a plain solve.  A plain
# rxSolve() queries each row once, so rxode2's own suite cannot reach the
# fault; only a fit can.  rxode2 carries bench/lincmt_phi_fit_engage_check.R
# as the reproducer, and this is the durable regression test.
#
# Which matrix: rxode2's linCmtSensPhi selects between three orderings of the
# same closed-form solution.  `1` assembles the matrix by probing the kernel
# with unit-basis prior states -- one evaluation per direction per column, so
# it is only worth paying where the interval repeats, and it is the mode the
# engage rule above governs and the `phiBuild`/`phiRows` counters measure.
# `TRUE`/`2` (rxode2's DEFAULT since the closed-form assembly landed) costs
# about one kernel evaluation, so it needs no evidence of recurrence at all
# and is counted separately as `phiAnalyticRows`.  This file therefore asks
# for `1` explicitly: passing `TRUE` silently moved it onto the analytic path,
# where `phiBuild` is 0 by construction and the engage rule it exists to pin
# was not being exercised at all.
#
# Skips cleanly on an rxode2 without the mechanism.

# The model function is an rxode2 ini()/model() DSL block, which lintr's
# object_usage cannot follow (every assignment looks unused).
# nolint start: object_usage_linter.
.phiMod <- function() {
  ini({
    lka <- log(1.2)
    lcl <- log(4)
    lv <- log(30)
    eta.cl ~ 0.1
    eta.v ~ 0.1
    prop.sd <- 0.2
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl) * exp(eta.cl)
    v <- exp(lv) * exp(eta.v)
    cp <- linCmt()
    cp ~ prop(prop.sd)
  })
}
# nolint end

# Does the loaded rxode2 have the transition-matrix path and its counters?
.phiFitCapable <- function() {
  .ns <- asNamespace("rxode2")
  if (!exists("linCmtSeqStats", envir = .ns, inherits = FALSE)) {
    return(FALSE)
  }
  # rxControl() takes ..., so the argument is declared on rxSolve().
  if (!("linCmtSensPhi" %in% names(formals(rxode2::rxSolve)))) {
    return(FALSE)
  }
  "phiBuild" %in% names(utils::getFromNamespace("linCmtSeqStats", "rxode2")(FALSE))
}

.phiStats <- function(reset = FALSE) {
  utils::getFromNamespace("linCmtSeqStats", "rxode2")(reset)
}

# Deterministic observations; the fit only has to exercise the inner problem,
# so no simulation model is compiled here.
.phiDat <- function(times, nid = 4L) {
  do.call(rbind, lapply(seq_len(nid), function(i) {
    cp <- 100 / 30 * (exp(-0.13 * times) - exp(-1.2 * times)) * (1 + 0.05 * i)
    rbind(
      data.frame(id = i, time = 0, amt = 100, evid = 1, dv = NA_real_),
      data.frame(
        id = i, time = times, amt = 0, evid = 0,
        dv = cp * (1 + 0.1 * sin(seq_along(times)))
      )
    )
  }))
}

# One posthoc fit, returning the counters it accumulated and its objective.
# phi = 1L is the probe-assembled matrix the engage rule governs; see the
# header.  FALSE is the row-by-row ordering, TRUE/2 the analytic assembly.
.phiFit <- function(times, phi = 1L) {
  .dat <- .phiDat(times)
  invisible(.phiStats(TRUE))
  .rx <- rxode2::rxControl(
    cores = 1L, linCmtSensType = "AD",
    linCmtSensPhi = phi
  )
  .ctl <- nlmixr2est::foceiControl(
    maxOuterIterations = 0L, print = 0L,
    calcTables = FALSE, covMethod = "", rxControl = .rx
  )
  .f <- suppressMessages(nlmixr2est::nlmixr2(.phiMod, .dat, "focei", .ctl))
  list(stats = .phiStats(TRUE), objf = .f$objf)
}

# Gaps strictly increasing: no interval ever recurs.
.phiTimesNonRepeating <- cumsum(seq(0.30, 0.70, length.out = 24))
# One repeating gap: the regime the matrix exists for.
.phiTimesRepeating <- seq(1, 24, by = 1)

test_that("intervals that never repeat build no transition matrix in a fit", {
  skip_on_cran()
  skip_if_not(.phiFitCapable())
  .r <- .phiFit(.phiTimesNonRepeating)
  # The regression this pins: a re-query of the same row is not evidence of
  # recurrence, so nothing may be assembled here.
  expect_equal(.r$stats[["phiBuild"]], 0L)
  expect_equal(.r$stats[["phiRows"]], 0L)
})

test_that("a repeating design engages the matrix and reuses it across rows", {
  skip_on_cran()
  skip_if_not(.phiFitCapable())
  .r <- .phiFit(.phiTimesRepeating)
  expect_gt(.r$stats[["phiBuild"]], 0L)
  # Reuse measured at 92 rows per build; a conservative floor catches a
  # collapse back to per-row assembly without pinning the exact ratio.
  expect_gt(.r$stats[["phiRows"]], 10L * .r$stats[["phiBuild"]])
  # The DEFAULT is the analytic assembly, which needs no recurrence and so
  # reports no builds -- cover it here too, or a default that stopped
  # engaging entirely would show up nowhere.
  .a <- .phiFit(.phiTimesRepeating, phi = TRUE)
  expect_equal(.a$stats[["phiBuild"]], 0L)
  expect_gt(.a$stats[["phiAnalyticRows"]], 0L)
  expect_equal(.a$objf, .r$objf, tolerance = 1e-8)
})

test_that("engaging the transition matrix does not move the objective", {
  skip_on_cran()
  skip_if_not(.phiFitCapable())
  .on <- .phiFit(.phiTimesRepeating, phi = 1L)
  .off <- .phiFit(.phiTimesRepeating, phi = FALSE)
  # linCmtSensPhi=FALSE must leave the matrix path entirely unused.
  expect_equal(.off$stats[["phiBuild"]], 0L)
  expect_gt(.on$stats[["phiBuild"]], 0L)
  # The two are the same exact closed-form solution evaluated in different
  # operation orders, so they agree to floating-point round-off.
  expect_equal(.on$objf, .off$objf, tolerance = 1e-8)
})
