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

  test_that("foce=\"foce+\" equals \"nonmem\" for pure additive error", {
    fn <- suppressWarnings(suppressMessages(
      nlmixr(add.cmt, d, "foce",
             foceiControl(print = 0L, calcTables = FALSE, outerOpt = "nlminb", covMethod = ""))))
    fp <- suppressWarnings(suppressMessages(
      nlmixr(add.cmt, d, "foce",
             foceiControl(print = 0L, calcTables = FALSE, outerOpt = "nlminb", covMethod = "", foce = "foce+"))))
    expect_equal(fn$objective, fp$objective, tolerance = 1e-4)
  })

  test_that("focep/mfocep/ifocep equal foce with foce=\"foce+\"", {
    ## foce+ reference for the focep alias check, cached as a value: focep IS
    ## foce + foce="foce+", so this arm is the entry point being aliased, not the
    ## behaviour under test, and foce+ itself is exercised in the two tests above.
    ## See helper-gradref.R.
    ref <- .numRef("fit-focep-ref-focepluss", function()
      list(objective = suppressWarnings(suppressMessages(
        nlmixr(one.cmt, d, "foce",
               foceiControl(print = 0L, calcTables = FALSE, outerOpt = "nlminb",
                            covMethod = "", foce = "foce+"))))$objective))
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
    # Measured under the CURRENT defaults (foce="nonmem" default, resetThetaP=0):
    #   foce+ / focep      122.682151  (identical -- focep IS foce + foce="foce+")
    #   mfocep             118.243871
    #   ifocep             118.202575
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
    expect_equal(fM$objective, 118.243871, tolerance = 1e-3)
    expect_equal(fI$objective, 118.202575, tolerance = 1e-3)
    # ... and keep them in the same neighbourhood as the plain fit.
    expect_lt(abs(fM$objective - ref$objective), 6)
    expect_lt(abs(fI$objective - ref$objective), 6)
  })
})
