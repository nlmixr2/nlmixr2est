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

  test_that("a flat random effect contributes nothing to the objective", {
    # A mu-referenced parameter whose omega is declared zero carries no
    # between-subject variability, so it must contribute no `eta^2/omega` to
    # the quadratic form and no `omega` to the log-determinant -- it is taken
    # out of the Cholesky (foceiOmegaDropFlat(), src/inner.cpp).  The variance
    # left in the omega matrix for it is a placeholder that keeps the matrix
    # invertible, and the test of "it never reaches the likelihood" is that the
    # objective does not depend on it.
    #
    # Verified directly while developing this by rebuilding with the internal
    # placeholder at 1 and at 100: focei 290.6220260211 and impmap
    # 193.5933442926 both times, to every digit.  What can be asserted from
    # here is the same invariance through the public interface -- the declared
    # value varies, the objective does not.
    .mk <- function(v) {
      .txt <- sprintf('function() {
        ini({ tka <- 0.45; tcl <- 1; tv <- 3.45
              eta.mu.tka ~ fix(%s); eta.mu.tv ~ fix(%s)
              eta.iiv ~ 0.1; add.sd <- 0.7 })
        model({ ka <- exp(tka + eta.mu.tka); cl <- exp(tcl + eta.iiv)
                v <- exp(tv + eta.mu.tv); linCmt() ~ add(add.sd) }) }', v, v)
      eval(parse(text = .txt))
    }
    .ofv <- function(m, est) {
      .ctl <- if (est == "focei") {
        foceiControl(print = 0, covMethod = "", maxOuterIterations = 0)
      } else {
        impmapControl(nIter = 3L, isample = 50L, print = 0, covMethod = "")
      }
      suppressWarnings(nlmixr2(m, nlmixr2data::theo_sd, est = est,
                               control = .ctl))$objf
    }
    for (.est in c("focei", "impmap")) {
      expect_equal(.ofv(.mk("0"), .est), .ofv(.mk("1e-9"), .est),
                   tolerance = 1e-10, info = .est)
    }
  })

  test_that("the exploration value is connected, and stays out of the omega", {
    # saemControl(zeroOmegaTune=) sets how far a flat random effect explores.
    #
    # The obvious assertion -- that the objective does not depend on it -- is
    # the WRONG test, and it passed for a long time while the parameter was
    # disconnected entirely: the pre-processing hook rewrote the omega before
    # the control was attached, so every fit used the default and 0.1 vs 4 gave
    # bit-identical results to every digit.  An inert input produces exactly
    # the same green test as a held invariant.
    #
    # So what is asserted is that the value REACHES the sampler (it changes the
    # trajectory) while remaining absent from what the model claims (omega is
    # reported as the declared zero).  Whether the declared zero contributes to
    # the objective is a separate property, tested above via fix(0)/fix(1e-9).
    .m <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.mu.tka ~ fix(0)
        eta.mu.tv ~ fix(0)
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
    .fit <- function(tune) {
      suppressWarnings(nlmixr2(.m(), nlmixr2data::theo_sd, est = "saem",
                               control = saemControl(nBurn = 30, nEm = 30,
                                                     print = 0, covMethod = "",
                                                     zeroOmegaTune = tune)))
    }
    .a <- .fit(0.1)
    .b <- .fit(4)
    # connected: a 40x change in sampling width cannot leave the fit identical
    expect_false(isTRUE(all.equal(.a$objf, .b$objf, tolerance = 1e-12)))
    # and still not a parameter of the model, at either width
    for (.f in list(.a, .b)) {
      expect_equal(.f$omega["eta.mu.tka", "eta.mu.tka"], 0)
      expect_equal(.f$omega["eta.mu.tv", "eta.mu.tv"], 0)
    }
  })
  test_that("the two opt-in M-step options validate and round-trip", {
    .c <- saemControl(zeroOmegaAnneal = 1000, zeroOmegaDirect = TRUE)
    expect_equal(.c$zeroOmegaAnneal, 1000)
    expect_true(.c$zeroOmegaDirect)
    # both default OFF -- they change the M-step, so they are opt-in
    .d <- saemControl()
    expect_equal(.d$zeroOmegaAnneal, 0)
    expect_false(.d$zeroOmegaDirect)
    expect_error(saemControl(zeroOmegaAnneal = -1), "zeroOmegaAnneal")
    expect_error(saemControl(zeroOmegaDirect = 1), "zeroOmegaDirect")
    .i <- impmapControl(zeroOmegaDirect = TRUE, zeroOmegaMaxEval = 10L)
    expect_true(.i$zeroOmegaDirect)
    expect_equal(.i$zeroOmegaMaxEval, 10L)
    expect_false(impmapControl()$zeroOmegaDirect)
    expect_error(impmapControl(zeroOmegaMaxEval = 0), "zeroOmegaMaxEval")
  })

  test_that("saem still moves the theta under each M-step option", {
    # The claim these options exist to fix is that the omega^-1-weighted GLS
    # cannot move such a theta at all; whatever route is taken, the theta must
    # not come back as its ini() value, and the declared zero must survive.
    .run <- function(...) {
      suppressWarnings(nlmixr2(.helperMod(), nlmixr2data::theo_sd, est = "saem",
                               control = saemControl(nBurn = 30, nEm = 30,
                                                     print = 0, covMethod = "",
                                                     ...)))
    }
    for (.nm in c("anneal", "direct", "both")) {
      .f <- switch(.nm,
                   anneal = .run(zeroOmegaAnneal = 1000),
                   direct = .run(zeroOmegaDirect = TRUE),
                   both   = .run(zeroOmegaAnneal = 1000, zeroOmegaDirect = TRUE))
      .est <- unname(.f$parFixedDf["tka", "Estimate"])
      expect_true(is.finite(.est), info = .nm)
      expect_false(isTRUE(all.equal(.est, 0.45, tolerance = 1e-7)), info = .nm)
      expect_equal(.f$omega["eta.mu.tka", "eta.mu.tka"], 0, info = .nm)
      expect_true(.f$omega["eta.iiv", "eta.iiv"] > 0, info = .nm)
    }
  })

  test_that("imp already moves a flat-omega mu theta, with or without the option", {
    # NOT the saem defect.  impMuInterceptStep() is theta += mean(eta), which
    # looks like the same degenerate update -- but imp gets its etas from the
    # inner MAP, and a flat random effect is dropped from Omega^-1 entirely
    # (foceiOmegaDropFlat()), so nothing penalizes it and the MAP moves it
    # freely against the data.  mean(eta) is therefore informative and the
    # theta does move.  saem cannot do this because it SAMPLES the eta from a
    # prior whose variance is ~0, so the column never leaves its prior mean.
    #
    # zeroOmegaDirect is still offered for imp (same NONMEM eqs. 1.47-1.52
    # route), but it is a refinement here rather than a repair -- measured on
    # this fixture it agreed with the mean-shift to every digit.
    .run <- function(...) {
      suppressWarnings(nlmixr2(.helperMod(), nlmixr2data::theo_sd, est = "impmap",
                               control = impmapControl(nIter = 5L, isample = 50L,
                                                       print = 0, covMethod = "",
                                                       ...)))
    }
    .off <- .run()
    .on <- .run(zeroOmegaDirect = TRUE, zeroOmegaMaxEval = 15L)
    for (.f in list(off = .off, on = .on)) {
      .est <- unname(.f$parFixedDf["tka", "Estimate"])
      expect_true(is.finite(.est))
      expect_false(isTRUE(all.equal(.est, 0.45, tolerance = 1e-7)))
      # the declared zero still reaches the report either way
      expect_equal(.f$omega["eta.mu.tka", "eta.mu.tka"], 0)
    }
  })

})
