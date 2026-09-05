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

  test_that("the exploration value stays out of the objective", {
    # saemControl(zeroOmegaTune=) sets how far a flat random effect explores so
    # the M-step has a conditional mean to shift its theta by.  It is machinery,
    # not a model parameter, so it must not reach the likelihood -- which is the
    # same defect a working variance in Omega would have.
    #
    # `fix(0)` vs `fix(1e-9)` cannot detect this: both spellings are rewritten
    # to the same tuning value for saem, so they agree either way.  Varying the
    # tuning value itself is what tests it.
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
    .ofv <- function(tune) {
      suppressWarnings(nlmixr2(.m(), nlmixr2data::theo_sd, est = "saem",
                               control = saemControl(nBurn = 30, nEm = 30,
                                                     print = 0, covMethod = "",
                                                     zeroOmegaTune = tune)))$objf
    }
    expect_equal(.ofv(0.1), .ofv(0.5), tolerance = 1e-10)
  })
})
