test_that("focei completes at cores>=2 without heap corruption (Windows cross-DLL OpenMP)", {
  # Regression test for the Windows-only heap corruption that crashed
  # windows-latest R-CMD-check while building focei fits.
  #
  # FOCEi's inner loop (innerOpt(), inner.cpp) runs `#pragma omp parallel
  # for` and each worker thread calls into rxode2's per-subject
  # ind_solve().  On Windows/Rtools rxode2 and nlmixr2est each statically
  # link their OWN libgomp, so inside rxode2 omp_get_thread_num() returns 0
  # for every one of nlmixr2est's worker threads.  All of rxode2's
  # per-thread solve buffers (the LSODA context pool, glhs/gon/solveSave
  # scratch, tolerance and dose arrays, linCmt scratch) then collapse onto
  # slot 0 and concurrent workers write / lsoda_free the same memory ->
  # heap corruption -> segfault.
  #
  # The fix hands rxode2 the real thread id via rxode2::setRxThreadId()
  # (called around each per-subject solve in innerOpt) so the per-thread
  # indexing is correct.
  #
  # The corruption is benign on Linux/macOS (one shared libgomp -> the
  # thread ids are already consistent) and at cores==1, so it only fires on
  # Windows at cores>=2.  It must run in a fresh callr subprocess at
  # cores>=2: if the corruption returns the child segfaults and callr::r()
  # raises an error -> test failure.  maxOuterIterations=0L keeps it cheap
  # while still exercising the full parallel inner solve.
  skip_on_cran()
  skip_if_not_installed("callr")

  obj <- callr::r(
    function() {
      library(nlmixr2est)
      rxode2::setRxThreads(2L)
      one.compartment <- function() {
        ini({
          tka <- log(1)
          tcl <- log(10)
          tv <- log(35)
          eta.ka ~ 0.1
          eta.cl ~ 0.1
          eta.v ~ 0.1
          add.sd <- 0.1
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
      fit <- nlmixr2(
        one.compartment, nlmixr2data::theo_sd, est = "focei",
        control = foceiControl(print = 0L, maxOuterIterations = 0L)
      )
      as.numeric(fit$objDf$OBJF[1])
    }
  )
  # If the fit completed the OBJ is finite; pin the value (computed at the
  # initial estimates, maxOuterIterations=0L) so a non-crashing but raced
  # result is also caught.
  expect_true(is.finite(obj))
  expect_equal(obj, 4613.033, tolerance = 1e-3)
})
