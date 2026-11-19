# Guards the bit-identical, thread-count-independent parallelization of the imp
# family: the E-step, the E-step proposal build, the theta-score M-step, and the
# Monte-Carlo covariance are all parallelized over subjects, so a full fit must be
# byte-for-byte identical at any thread count.  Run on every CI runner,
# Windows' CRAN-only subset included, where a stale rxode2 thread id made the
# linCmt() fits depend on the thread count; not on CRAN itself.
nmTest(
  {
    test_that("imp family full fit is thread-count independent (E/M-step + cov)", {
      .skipOnCranOnly()
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
      # covMethod="imp" (default) exercises the parallel Monte-Carlo covariance;
      # nIter>=3 exercises several E/M-step iterations.  imp (no MAP search) and
      # impmap (MAP proposal build) cover both proposal paths.
      .run <- function(est, nthr) {
        rxode2::setRxThreads(nthr)
        rxode2::rxSetSeed(42)
        .ctl <- if (est == "imp") {
          impControl(print = 0L, nIter = 3L, isample = 200L)
        } else {
          impmapControl(print = 0L, nIter = 3L, isample = 200L)
        }
        .f <- suppressWarnings(nlmixr2(one.cmt, .dat, est, .ctl))
        list(objf = .f$objf, fixef = .f$fixef, omega = .f$omega, cov = .f$cov)
      }
      for (.est in c("imp", "impmap")) {
        .r1 <- .run(.est, 1L)
        .r4 <- .run(.est, 4L)
        expect_identical(.r1$objf, .r4$objf)
        expect_identical(.r1$fixef, .r4$fixef)
        expect_identical(.r1$omega, .r4$omega)
        expect_identical(.r1$cov, .r4$cov)
      }
    })
  },
  skipOnCran = FALSE
)
