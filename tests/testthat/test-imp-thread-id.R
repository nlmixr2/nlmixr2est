# An ODE model on three subjects: with two threads, each thread solves a
# second subject, the case where a stale rxode2 thread id put two threads on
# one liblsoda slot (a crash on Windows, where rxode2's own
# omp_get_thread_num() is 0 for every nlmixr2est thread)
.impThreadMod <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.cl ~ 0.3
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

.impThreadData <- function() {
  .d <- nlmixr2data::theo_sd
  .d[.d$ID %in% 1:3, ]
}

# c(checked, miss, openmp): the parallel subject iterations checked, and those
# that began on another thread's rxode2 id, since the last reset; and whether
# nlmixr2est has OpenMP
.impThreadMiss <- function(reset = FALSE) {
  .Call(`_nlmixr2est_nmRxThreadIdMissTest`, reset, PACKAGE = "nlmixr2est")
}

# Turn the thread-id check on or off; returns whether it was on
.impThreadTrack <- function(on) {
  .Call(`_nlmixr2est_nmRxThreadIdTrack`, on, PACKAGE = "nlmixr2est")
}

.impThreadFit <- function(est, nthr) {
  rxode2::setRxThreads(nthr)
  suppressMessages(suppressWarnings(
    nlmixr2(.impThreadMod, .impThreadData(), est, control = list(print = 0L))
  ))
}

# Run on every CI runner, Windows' CRAN-only subset included, where the stale
# id crashed; not on CRAN itself
nmTest(
  {
    test_that("imp and impmap give each thread its own rxode2 id for every subject", {
      .skipOnCranOnly()
      .old <- rxode2::getRxThreads()
      withr::defer(rxode2::setRxThreads(.old))
      # setRxThreads() clamps to the CPUs and OMP_THREAD_LIMIT; with one
      # thread the region does not fork and nothing is checked
      rxode2::setRxThreads(2L)
      skip_if(rxode2::getRxThreads() < 2L, "one thread")
      skip_if(.impThreadMiss()[["openmp"]] == 0L, "nlmixr2est has no OpenMP")
      .wasOn <- .impThreadTrack(TRUE)
      withr::defer(.impThreadTrack(.wasOn))
      .want <- c(imp = 18.5322024276, impmap = 18.5303802572)
      for (.est in names(.want)) {
        .one <- .impThreadFit(.est, 1L)
        .impThreadMiss(TRUE)
        .two <- .impThreadFit(.est, 2L)
        .m <- .impThreadMiss()
        # the region forked, so the check ran
        expect_gt(.m[["checked"]], 0L, label = .est)
        expect_identical(.m[["miss"]], 0L, label = .est)
        # the objective does not depend on the thread count
        expect_equal(.two$objective, .one$objective, tolerance = 1e-10, label = .est)
        expect_equal(.two$objective, .want[[.est]], tolerance = 1e-6, label = .est)
      }
    })
  },
  skipOnCran = FALSE
)
