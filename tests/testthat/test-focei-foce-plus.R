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
    ref <- suppressWarnings(suppressMessages(
      nlmixr(one.cmt, d, "foce",
             foceiControl(print = 0L, calcTables = FALSE, outerOpt = "nlminb", covMethod = "", foce = "foce+"))))
    fit <- suppressWarnings(suppressMessages(
      nlmixr(one.cmt, d, "focep",
             foceiControl(print = 0L, calcTables = FALSE, outerOpt = "nlminb", covMethod = ""))))
    expect_true(is.finite(fit$objective))
    expect_equal(fit$objective, ref$objective, tolerance = 1e-3)
    # The mu-profiled variants are NOT expected to match the plain fit here:
    # the foce+ per-subject inner problem is multi-modal, so the three land in
    # different conditional basins.  They must agree with EACH OTHER.
    #
    # This test used to additionally assert `fM <= ref + 0.5` -- that the
    # mu-group regression warm-starts the inner optimizer into DEEPER modes than
    # the plain fit reaches.  That held only while the plain fit was handicapped
    # by the inner eta reset discarding its converged EBEs.  With that fixed the
    # plain fit reaches the deepest mode on this fixture and the inequality
    # flips, so the directional claim is no longer a property of the estimators:
    #
    #                        objective   cold-start at its own estimates
    #   before the eta fix:
    #     foce+ / focep       120.9624   120.7041
    #     mfocep / ifocep     118.2396   118.3629
    #   after:
    #     foce+ / focep       116.6300   116.2202   <- genuinely deeper
    #     mfocep / ifocep     118.2396   118.3629   <- bit-identical; the mu
    #                                                  family is untouched by
    #                                                  the eta-reset fix
    #
    # So the mu variants are now the ones sitting in the shallower basin, which
    # is worth a look on its own (they take a different inner path and did not
    # benefit from the fix) but is not something this test should encode as an
    # expectation.  warm="save" (self-init inner Hessian) is still pinned: the
    # default warm="calc" recalculates the eta Hessian at the mu-regression's
    # restarted theta/eta and moves this fixture again.
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
    expect_equal(fM$objective, 118.2396, tolerance = 1e-3)
    expect_equal(fI$objective, 118.2399, tolerance = 1e-3)
    # ... and keep them in the same neighbourhood as the plain fit.
    expect_lt(abs(fM$objective - ref$objective), 3)
    expect_lt(abs(fI$objective - ref$objective), 3)
  })
})
