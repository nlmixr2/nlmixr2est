nmTest({
  test_that("trustControl() defaults and validation", {
    expect_equal(trustControl()$fterm, 1e-3)
    expect_equal(trustControl()$mterm, 1e-3)
    expect_equal(trustControl(sigdig = 4)$fterm, 1e-4)
    expect_equal(trustControl(fterm = 0.5)$mterm, 0.5)
    expect_equal(trustControl(fterm = 0.5, mterm = 0.1)$mterm, 0.1)
    expect_null(trustControl()$rinit)
    expect_null(trustControl()$rmax)
    expect_equal(trustControl(rmax = 2)$rmax, 2)

    expect_error(trustControl(rinit = 0))
    expect_error(trustControl(rmax = 0))
    expect_error(trustControl(rinit = 5, rmax = 1))

    .ctl <- trustControl(rinit = 0.1, rmax = 1)
    expect_equal(do.call(trustControl, .ctl)$rinit, 0.1)
    expect_equal(do.call(trustControl, .ctl)$rmax, 1)
  })

  test_that("trust registers in the nlm-family method listing", {
    .types <- nlmixr2AllEstType()
    expect_true("trust" %in% .types$est)
    expect_equal(.types$type[.types$est == "trust"], "Optimizer (NLM family)")
  })

  .oneCmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(centr) <- ka * depot - cl / v * centr
      cp <- centr / v
      cp ~ add(add.sd)
    })
  }

  test_that("est='trust' converges close to bobyqa", {
    skip_on_cran()
    .fB <- .nlmixr(.oneCmt, nlmixr2data::theo_sd, est = "bobyqa",
                   control = bobyqaControl(print = 0L, calcTables = FALSE))
    .fT <- .nlmixr(.oneCmt, nlmixr2data::theo_sd, est = "trust",
                   control = trustControl(print = 0L, calcTables = FALSE))
    .nT <- .nTrustOuter()

    expect_true(is.finite(.fB$objective))
    expect_true(is.finite(.fT$objective))
    expect_equal(.fT$objective, .fB$objective, tolerance = 1e-2)
    expect_equal(unname(.fT$theta), unname(.fB$theta), tolerance = 1e-2)

    # Positive evidence the C++-resident loop actually ran -- not just that
    # the numbers happen to agree (mirrors test-focei-trust-inner.R's own
    # .nTrustInner() convention for the same reason).
    expect_true(.nT > 0L)
  })

  test_that("est='trust' returnTrust=TRUE gives the raw trust_solve_c() output", {
    skip_on_cran()
    .raw <- .nlmixr(.oneCmt, nlmixr2data::theo_sd, est = "trust",
                    control = trustControl(print = 0L, calcTables = FALSE,
                                           returnTrust = TRUE))
    expect_true(is.list(.raw))
    expect_true(is.logical(.raw$converged))
    expect_true(is.numeric(.raw$value))
    expect_equal(dim(.raw$hessian), c(4L, 4L))
  })

  test_that("est='trust' reuses its own Hessian, skipping nlmixr2est's post-fit FD Hessian", {
    skip_on_cran()
    .calledTrust <- FALSE
    testthat::with_mocked_bindings(
      nlmixr2Hess = function(...) {
        .calledTrust <<- TRUE
        stop("nlmixr2Hess should not be called for est='trust'")
      },
      .package = "nlmixr2est",
      {
        .fT <- .nlmixr(.oneCmt, nlmixr2data::theo_sd, est = "trust",
                       control = trustControl(print = 0L, calcTables = FALSE))
      }
    )
    expect_false(.calledTrust)
    expect_true(is.finite(.fT$objective))

    # Control: the SAME mock trips for bobyqa (which has no Hessian of its
    # own and needs the post-fit FD recompute), proving the mock itself is
    # capable of detecting the call rather than being coincidentally unused.
    .calledBobyqa <- FALSE
    testthat::with_mocked_bindings(
      nlmixr2Hess = function(...) {
        .calledBobyqa <<- TRUE
        stop("nlmixr2Hess should be called for est='bobyqa'")
      },
      .package = "nlmixr2est",
      {
        expect_error(
          .nlmixr(.oneCmt, nlmixr2data::theo_sd, est = "bobyqa",
                  control = bobyqaControl(print = 0L, calcTables = FALSE))
        )
      }
    )
    expect_true(.calledBobyqa)
  })

  test_that("est='trust' warns when the Newton decrement contradicts trust_solve_c()'s own converged flag", {
    # The false-positive-convergence case this warning guards against (RcppTrust
    # reports converged==TRUE off its own internal step-size tolerance while the
    # true Newton step at that point is still large, verified in nlmTrustFit(),
    # src/nlm.cpp) is a genuine numerical edge case, not reliably reproducible
    # on demand from a real fit -- test the warning logic directly instead.
    expect_warning(
      nlmixr2est:::.trustWarnUnderConverged(list(underConverged = TRUE)),
      "Newton step|verified stationary point"
    )
    expect_silent(nlmixr2est:::.trustWarnUnderConverged(list(underConverged = FALSE)))
    expect_silent(nlmixr2est:::.trustWarnUnderConverged(list()))
  })
})
