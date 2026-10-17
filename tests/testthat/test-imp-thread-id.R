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

# The parallel subject iterations that began on another thread's rxode2 id
.impThreadMiss <- function(reset = FALSE) {
  .Call(`_nlmixr2est_nmRxThreadIdMissTest`, reset, PACKAGE = "nlmixr2est")
}

nmTest({
  test_that("imp and impmap give each thread its own rxode2 id for every subject", {
    skip_on_cran()
    .old <- rxode2::getRxThreads()
    rxode2::setRxThreads(2L)
    withr::defer(rxode2::setRxThreads(.old))
    # the objective does not depend on the thread count (1, 2 and 3 threads
    # give these to 10 digits)
    .want <- c(imp = 18.5322024276, impmap = 18.5303802572)
    for (.est in names(.want)) {
      .impThreadMiss(TRUE)
      .fit <- suppressMessages(suppressWarnings(
        nlmixr2(.impThreadMod, .impThreadData(), .est, control = list(print = 0L))
      ))
      expect_equal(.fit$objective, .want[[.est]], tolerance = 1e-9, label = .est)
      expect_identical(.impThreadMiss(), 0L, label = .est)
    }
  })
})
