test_that("nlmSetup does not segfault on first estimator call (no random effects)", {
  # Regression test for nlm.cpp:146 where
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
  # Pre-fix, this call would crash the test R process.  Post-fix the
  # call completes cleanly.  We use stats::nlm so the test does not
  # need minqa, and keep iterlim minimal so this stays cheap — the
  # crash is in setup, not in the actual optimization, so a single
  # iteration is enough to exercise the regression path.
  skip_on_cran()

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
  expect_s3_class(fit, "nlmixr2FitCore")
})
