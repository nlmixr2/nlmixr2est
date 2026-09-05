nmTest({
  # "$OMEGA (0.0 FIXED)" on a MU_ helper eta is how a NONMEM control stream
  # mu-references a plain theta: the random effect carries no variability of its
  # own, it exists so the EM has a conditional mean to shift the theta by.
  #
  # Held the way an ordinary fix() is held, that freezes everything -- saem's
  # mu-theta M-step is weighted by omega^-1, so the column gets infinite weight
  # and cannot move off its prior mean, which IS the current theta.  Measured on
  # the model below before this was handled: tka came back as 0.4500, exactly
  # its ini() value, reported as an estimate.
  #
  # NB one zero-variance helper on its own does not reproduce it -- the near-PD
  # correction on Gamma2_phi1 rescues a single zero diagonal.  It takes more
  # than one, which is also the shape a ported control stream has (one helper
  # per distribution parameter).
  .helperMod <- function() {
    .f <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.mu.tka ~ fix(1e-8)
        eta.mu.tv ~ fix(1e-8)
        eta.iiv ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.mu.tka)
        cl <- exp(tcl + eta.iiv)
        v <- exp(tv + eta.mu.tv)
        linCmt() ~ add(add.sd)
      })
    }
    .f()
  }

  test_that("a mu-referenced variance declared as zero is detected", {
    expect_setequal(nlmixr2est:::.zeroOmegaMuRefEtas(.helperMod()),
                    c("eta.mu.tka", "eta.mu.tv"))
  })

  test_that("it does not fire on an ordinary fix()ed variance", {
    .m <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ fix(0.3)
        eta.cl ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    expect_equal(length(nlmixr2est:::.zeroOmegaMuRefEtas(.m())), 0L)
  })

  test_that("saem moves the theta and still reports the declared zero", {
    .f <- suppressWarnings(nlmixr2(.helperMod(), nlmixr2data::theo_sd,
                                   est = "saem",
                                   control = saemControl(nBurn = 30, nEm = 30,
                                                         print = 0,
                                                         covMethod = "")))
    # Moved: the whole point.  Frozen, this was 0.45 to every digit.
    expect_false(isTRUE(all.equal(unname(.f$parFixedDf["tka", "Estimate"]),
                                  0.45, tolerance = 1e-7)))
    expect_true(is.finite(.f$parFixedDf["tka", "Estimate"]))
    # And the working variance does not escape into the report -- it was
    # machinery, and reporting it would claim between-subject variability the
    # model does not have.  (This half failed silently at first: the marker was
    # stashed with rxAssignControlValue(), which does not survive to the fit.)
    for (.e in c("eta.mu.tka", "eta.mu.tv")) {
      expect_equal(.f$omega[.e, .e], 0, info = .e)
      expect_true(all(.f$omega[.e, ] == 0), info = .e)
    }
    # The real random effect is untouched.
    expect_true(.f$omega["eta.iiv", "eta.iiv"] > 0)
  })
})
