nmTest({
  # #1073: the reported omega was snapshotted into Gamma2_phi1Report BEFORE the
  # fix()ed-value restore (and before the variance floors), so a fix()ed eta
  # variance came back as the M-step's unconstrained estimate.  The sampler was
  # always right -- only the report was wrong, which is the hard kind to notice.

  .fixMuMod <- function() {
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

  # the same thing for an eta that is NOT `theta + eta`, so it is carried
  # through nonMuEtas rather than mu-referenced -- the drift was never about
  # how the eta is parameterized
  .fixNonMuMod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      rxz.eta.ka ~ fix(1)
      eta.cl ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + 0.3 * logit(pnorm(rxz.eta.ka)))
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  .ctl <- function() {
    saemControl(nBurn = 20, nEm = 20, print = 0, seed = 42, calcTables = FALSE,
                covMethod = "")
  }

  test_that("a fix()ed eta variance is reported as the value it was fixed at", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.fixMuMod, nlmixr2data::theo_sd, "saem", .ctl())))
    expect_equal(unname(.fit$omega["eta.ka", "eta.ka"]), 0.3)
    # the estimated one still moved, so the fit is not simply echoing ini()
    expect_false(isTRUE(all.equal(unname(.fit$omega["eta.cl", "eta.cl"]), 0.1)))

    # The reported matrix and the one the sampler used must agree ON THE FIXED
    # CELL.  That comparison, not fit$omega alone, is what catches a
    # regression: fit$omega cannot tell "fixed at 0.3" from "estimated and
    # landed near 0.3".  Only the fixed cells are compared because the report
    # is deliberately left alone elsewhere -- a floored or PD-repaired variance
    # still reports its pre-repair value.
    .w <- which(diag(.fit$saem$Gamma2_phi1) == 0.3)
    expect_equal(length(.w), 1L)
    expect_equal(.fit$saem$Gamma2_phi1Report[.w, .w], .fit$saem$Gamma2_phi1[.w, .w])
  })

  test_that("the SA covariance phase does not undo the fixed report", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    # covMethod="sa" is the default, and it runs nSaCov extra iterations after
    # the fit; covMethod="" skips them entirely, so the other tests here never
    # see that phase at all.
    #
    # This does NOT exercise the _savGamma2_phi1Report snapshot/restore around
    # the phase: the gain is frozen at zero there (`pas`/`pash` are zero-padded,
    # src/saem.cpp), so the sufficient statistics and hence Gamma2_phi1 do not
    # move, and neutering the restore leaves both assertions below passing
    # (measured).  What it does establish is the user-visible invariant -- the
    # phase must not change the reported omega, fixed cell or otherwise.
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.fixMuMod, nlmixr2data::theo_sd, "saem",
              saemControl(nBurn = 20, nEm = 20, print = 0, seed = 42,
                          calcTables = FALSE, covMethod = "sa", nSaCov = 20L))))
    expect_equal(unname(.fit$omega["eta.ka", "eta.ka"]), 0.3)
    # ... and the estimated entries are untouched by the phase too, not just
    # the fixed one.
    .noCov <- suppressWarnings(suppressMessages(
      nlmixr2(.fixMuMod, nlmixr2data::theo_sd, "saem",
              saemControl(nBurn = 20, nEm = 20, print = 0, seed = 42,
                          calcTables = FALSE, covMethod = ""))))
    expect_equal(unname(.fit$omega), unname(.noCov$omega))
  })

  test_that("a non-mu-referenced eta's fix()ed variance is reported too", {
    skip_on_cran()
    skip_if_not_installed("nlmixr2data")
    .fit <- suppressWarnings(suppressMessages(
      nlmixr2(.fixNonMuMod, nlmixr2data::theo_sd, "saem", .ctl())))
    expect_equal(unname(.fit$omega["rxz.eta.ka", "rxz.eta.ka"]), 1)
    .w <- which(diag(.fit$saem$Gamma2_phi1) == 1)
    expect_equal(length(.w), 1L)
    expect_equal(.fit$saem$Gamma2_phi1Report[.w, .w], .fit$saem$Gamma2_phi1[.w, .w])
  })
})
