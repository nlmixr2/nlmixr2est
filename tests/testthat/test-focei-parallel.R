nmTest({
  one.compartment <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v  <- exp(tv  + eta.v)
      d/dt(depot)  <- -ka * depot
      d/dt(center) <-  ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  test_that("FOCEi cores=2 matches cores=1 baseline", {
    fit1 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "",
                                          rxControl = rxode2::rxControl(cores = 1L)))
    fit2 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "",
                                          rxControl = rxode2::rxControl(cores = 2L)))
    expect_s3_class(fit1, "nlmixr2FitData")
    expect_s3_class(fit2, "nlmixr2FitData")
    expect_equal(fit1$objective, fit2$objective, tolerance = 1e-4)
    expect_equal(as.numeric(fit1$theta), as.numeric(fit2$theta), tolerance = 1e-4)
  })

  test_that("FOCEi cores=2 matches cores=1 with mceta=0 (always reset)", {
    # mceta=0 path in innerOpt1: zero ETAs each call, no R API.
    fit1 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "", mceta = 0L,
                                          rxControl = rxode2::rxControl(cores = 1L)))
    fit2 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "", mceta = 0L,
                                          rxControl = rxode2::rxControl(cores = 2L)))
    expect_equal(fit1$objective, fit2$objective, tolerance = 1e-4)
  })

  test_that("FOCEi cores=2 with mceta>=1 (pre-drawn samples) does not crash", {
    # mceta >= 1 used to force cores=1 because the sampling path called the
    # R API.  Samples are now pre-drawn serially before the parallel region.
    # The objective need not match cores=1 bit-for-bit because the optimizer
    # trajectory differs across cores; only require that both fits complete.
    fit1 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "", mceta = 2L,
                                          rxControl = rxode2::rxControl(cores = 1L)))
    fit2 <- .nlmixr(one.compartment, nlmixr2data::theo_sd, est = "focei",
                    control = foceiControl(print = 0, covMethod = "", mceta = 2L,
                                          rxControl = rxode2::rxControl(cores = 2L)))
    expect_s3_class(fit1, "nlmixr2FitData")
    expect_s3_class(fit2, "nlmixr2FitData")
    expect_true(is.finite(fit1$objective))
    expect_true(is.finite(fit2$objective))
  })

  # NOTE: Models whose f() / dur() / rate() depend on ETAs (or any model
  # that triggers the shi21 finite-difference path) currently crash under
  # cores >= 2.  The crash is rooted in shared-global op->neq mutation in
  # shi21EtaGeneral racing with concurrent inner solves.  See the reprex
  # at inst/reprex_fbio_eta_parallel.R.  Fix path: rxode2 par_solve.cpp
  # needs to consult the per-individual neqOverride C-API
  # (already exposed at slots 56/57 in rxode2 5.0.x).  Tracking as a
  # follow-up issue; do not enable cores > 1 for those models until then.
})
