test_that("nlmSetup does not segfault on first estimator call in a fresh R session", {
  # Regression test for nlm.cpp where
  #   nlmOp.stickyRecalcN2Per.assign((size_t)getRxNsub(rx), 0)
  # was called BEFORE rxode2::rxSolve_() initializes the global `rx`
  # pointer.  In a fresh R session `rx` is still NULL and the read
  # inside getRxNsub() dereferences a NULL pointer ("memory not mapped,
  # address 0x10").
  #
  # Affected every pooled estimator that goes through .nlmSetupEnv:
  #   bobyqa, nlm, optim, nls, nlminb, lbfgsb3c, n1qn1, newuoa, uobyqa.
  # focei uses foceiSetup_ and was unaffected.
  #
  # The regression only fires when `rx` has not yet been populated by
  # an earlier rxode2 solve or focei fit, so this test MUST run in a
  # truly fresh R session -- running inline in the testthat process is
  # not enough because earlier tests will already have populated `rx`
  # (per Matt Fidler's review on #653).  callr spawns a clean
  # subprocess; if the segfault returns the child dies and
  # callr::r() raises a callr_status_error which surfaces as a test
  # failure.  stats::nlm keeps the test cheap; the crash is in setup
  # so iterlim = 1L exercises the regression path.
  skip_on_cran()
  skip_if_not_installed("callr")

  cls <- callr::r(
    function() {
      library(nlmixr2est)
      mod <- function() {
        ini({
          tka <- 0.45
          tcl <- log(2.7)
          tv <- 3.45
          add.sd <- 0.7
        })
        model({
          ka <- exp(tka)
          cl <- exp(tcl)
          v <- exp(tv)
          linCmt() ~ add(add.sd)
        })
      }
      fit <- nlmixr2(
        mod, nlmixr2data::theo_sd, est = "nlm",
        control = nlmControl(iterlim = 1L, print.level = 0L)
      )
      class(fit)
    }
  )
  expect_true("nlmixr2FitCore" %in% cls)
})
