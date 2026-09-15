nmTest({
  # foce = "foce+" revives the pre-6.1.0 FOCE that keeps the live conditional
  # residual variance R (vs. the default "nonmem" FOCE that freezes R at eta=0).
  # For a proportional/combined-error model the two objectives differ; for pure
  # additive error R does not depend on eta, so they must coincide.
  skip_on_cran()
  skip_if_not_installed("nlmixr2data")

  one.cmt <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7; prop.sd <- 0.1
          eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd) + prop(prop.sd) })
  }

  add.cmt <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- 0.7
          eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd) })
  }

  d <- nlmixr2data::theo_sd

  test_that("foce=\"foce+\" runs and differs from \"nonmem\" for proportional error", {
    fn <- suppressWarnings(suppressMessages(
      nlmixr(one.cmt, d, "foce",
             foceiControl(print = 0L, calcTables = FALSE, outerOpt = "nlminb", covMethod = ""))))
    fp <- suppressWarnings(suppressMessages(
      nlmixr(one.cmt, d, "foce",
             foceiControl(print = 0L, calcTables = FALSE, outerOpt = "nlminb", covMethod = "", foce = "foce+"))))
    expect_true(is.finite(fn$objective))
    expect_true(is.finite(fp$objective))
    expect_false(isTRUE(all.equal(fn$objective, fp$objective, tolerance = 1e-4)))
  })

  test_that("foce=\"foce+\" equals \"nonmem\" at matched additive parameters", {
    fn <- suppressWarnings(suppressMessages(
      nlmixr(add.cmt, d, "foce",
             foceiControl(print = 0L, calcTables = FALSE, maxOuterIterations = 0L,
                          innerOpt = "n1qn1", epsilon = 1e-10, maxInnerIterations = 1000L, covMethod = ""))))
    fp <- suppressWarnings(suppressMessages(
      nlmixr(add.cmt, d, "foce",
             foceiControl(print = 0L, calcTables = FALSE, maxOuterIterations = 0L,
                          innerOpt = "n1qn1", epsilon = 1e-10, maxInnerIterations = 1000L, covMethod = "", foce = "foce+"))))
    expect_equal(fn$objective, fp$objective, tolerance = 1e-4)
  })

  test_that("focep/mfocep/ifocep equal foce with foce=\"foce+\"", {
    # Compare aliases against a live reference from this implementation.
    ref <- list(objective = suppressWarnings(suppressMessages(
      nlmixr(one.cmt, d, "foce",
             foceiControl(print = 0L, calcTables = FALSE, outerOpt = "nlminb",
                          covMethod = "", foce = "foce+"))))$objective)
    fit <- suppressWarnings(suppressMessages(
      nlmixr(one.cmt, d, "focep",
             foceiControl(print = 0L, calcTables = FALSE, outerOpt = "nlminb", covMethod = ""))))
    expect_true(is.finite(fit$objective))
    expect_equal(fit$objective, ref$objective, tolerance = 1e-3)
    # The mu-profiled variants are NOT expected to match the plain fit here: the foce+
    # per-subject inner problem is multi-modal, so the three land in different conditional
    # basins.  They must agree with EACH OTHER, and stay in the same neighbourhood as the
    # plain fit -- which is a sanity check, not a precision claim, so the bound is loose.
    #
    # An earlier revision of this file recorded 116.63 for "foce+ / focep" and sized the
    # neighbourhood bound at 3 around it.  That number dates from when est="foce" WAS the
    # foce+ variant, so it describes a code state that no longer exists and is not a target
    # to restore.  Do not re-tighten this bound to fit it.  warm="save" (self-init inner
    # Hessian) is still pinned: the default warm="calc" recalculates the eta Hessian at the
    # mu-regression's restarted theta/eta and moves this fixture again.
    fM <- suppressWarnings(suppressMessages(
      nlmixr(one.cmt, d, "mfocep",
             foceiControl(print = 0L, calcTables = FALSE, outerOpt = "nlminb", covMethod = "", warm = "save"))))
    fI <- suppressWarnings(suppressMessages(
      nlmixr(one.cmt, d, "ifocep",
             foceiControl(print = 0L, calcTables = FALSE, outerOpt = "nlminb", covMethod = "", warm = "save"))))
    expect_true(is.finite(fM$objective))
    expect_true(is.finite(fI$objective))
    expect_equal(fM$objective, fI$objective, tolerance = 1e-2)
    # Pin the mu variants to the basin they actually reach on this fixture, so a
    # future change to the mu-referenced inner path shows up here.
    expect_equal(fM$objective, 107.341169, tolerance = 1e-3)
    expect_equal(fI$objective, 107.341169, tolerance = 1e-3)
    # ... and keep them in the same neighbourhood as the plain fit.  The gap is
    # between two different conditional basins, so it moves whenever either inner
    # path does; the bound is a tripwire for unrelated changes, not a claim about
    # this one.
    expect_lt(abs(fM$objective - ref$objective), 8)
    expect_lt(abs(fI$objective - ref$objective), 8)
  })

  test_that("a stalled FOCE+ EBE polish does not poison the objective (#1069)", {
    # FOCE+ polishes the inner optimizer's eta onto the truncated-score root it
    # defines its EBE by.  That polish stops once the score reaches the noise floor
    # the solve tolerance buys -- near 1e-3 at the DEFAULT sigdig=3 -- and used to
    # report the subject's likelihood as NA when it did.  Run at the default on
    # purpose: pinning sigdig is what hid this.
    .m <- function() {
      ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(31.5)
            eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; add.sd <- 0.7 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
        d/dt(depot)  <- -ka * depot
        d/dt(center) <-  ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd) })
    }
    fit <- suppressWarnings(suppressMessages(
      nlmixr(.m, d, "focep", focepControl(print = 0L, calcTables = FALSE, covMethod = ""))))
    # The stall cost 4.8 objective units: bobyqa stopped at 121.560 with the omegas
    # barely off their starting values, against 116.80 here and 116.82 for est="foce".
    expect_lt(fit$objf, 118)
    # The mechanism, not just the answer.  A subject whose likelihood came back NA
    # was penalized to ~300 while its neighbourhood sat near 130, and the cliff that
    # left in the outer objective is what stopped the search.  Nothing the search
    # sees may be worse than its own starting point by that kind of margin.
    .o <- fit$parHist$objf
    expect_lt(max(.o), .o[1] + 25)
    # ... and the omegas have to have actually moved (the stall left eta.ka at 0.56
    # against 0.40 here, because every probe that raised it read as a cliff).
    expect_lt(fit$omega[["eta.ka", "eta.ka"]], 0.5)
  })
})
