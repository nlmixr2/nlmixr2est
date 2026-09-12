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
