nmTest({
  test_that("imp's M-step reports its damping and mceta cannot blow up the thetas", {
    # The M-step used to take a raw Newton step on the non-mu structural thetas
    # (`solve(H, g)` guarded only by "solved" and "finite"), and a near-singular
    # IS-weighted Hessian satisfies both while returning an astronomically large
    # step.  One EM iteration then threw the whole parameter vector off the map,
    # and which mceta= the inner MAP used decided whether that happened.
    #
    # So the invariant is: changing only mceta= must not change the answer by
    # orders of magnitude.  The thetas here are deliberately non-mu-referenced
    # (ka enters through a non-linear covariate form) so the structural Newton
    # step is exercised rather than skipped.
    one <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.cl ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    ctl <- function(mceta) {
      impmapControl(nIter = 3L, isample = 50L, print = 0L, covMethod = "",
                    mceta = mceta)
    }
    f0 <- suppressWarnings(nlmixr2(one, nlmixr2data::theo_sd, est = "impmap",
                                   control = ctl(0L)))

    # The counters exist and a healthy fit needs no rescue.
    expect_true(is.numeric(f0$env$impMStepDamped))
    expect_true(is.numeric(f0$env$impMStepSkipped))
    expect_equal(f0$env$impMStepSkipped, 0)

    expect_true(is.finite(f0$objf))
    for (mc in c(10L, 100L)) {
      fm <- suppressWarnings(nlmixr2(one, nlmixr2data::theo_sd, est = "impmap",
                                     control = ctl(mc)))
      expect_true(is.finite(fm$objf))
      # Same model, same data, same iteration count -- only the inner MAP's
      # starting points differ, so the objective must stay in the same place.
      expect_true(abs(fm$objf - f0$objf) < 0.05 * abs(f0$objf))
      # And no theta may have walked off: exp(tcl) is a clearance of ~e.
      expect_true(all(is.finite(fm$parFixedDf$Estimate)))
      expect_true(all(abs(fm$parFixedDf$Estimate) < 50))
    }
  })
})

nmTest({
  test_that("an ill-conditioned M-step Hessian is damped, not solved approximately", {
    # The reproduction the damping exists for.  `tv` and `tv2` are both non-mu
    # structural thetas and enter only through a single sum, scaled so the pair
    # is numerically unidentifiable -- the M-step's theta-sensitivity Hessian is
    # then singular to working precision (measured rcond ~1e-21).
    #
    # Unguarded, `arma::solve()` answers that system with an approximate
    # least-squares solution AND prints
    #   "solve(): system is singular; rcond: 1.7998e-21; attempting approx solution"
    # to the console, once per M-step iteration.  What makes it a bug rather
    # than noise is that the approximate step was then taken: nothing in
    # "solved and finite" distinguishes it from a real Newton step.
    skip_on_cran()
    collinear <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        tv2 <- 0.01
        eta.cl ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + 0.001 * tv2)
        linCmt() ~ add(add.sd)
      })
    }
    f <- suppressWarnings(nlmixr2(collinear, nlmixr2data::theo_sd, est = "impmap",
                                  control = impmapControl(nIter = 5L, isample = 50L,
                                                          print = 0L, covMethod = "")))
    # The guard engaged on this model rather than sitting unused.
    expect_gt(f$env$impMStepDamped, 0)
    # ...and having engaged, it kept the thetas somewhere a model can live.
    expect_true(is.finite(f$objf))
    expect_true(all(is.finite(f$parFixedDf$Estimate)))
    expect_true(all(abs(f$parFixedDf$Estimate) < 50))
  })
})
