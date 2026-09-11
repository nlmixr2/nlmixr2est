nmTest({
  # #1073: the reporting snapshot `Gamma2_phi1Report` used to be taken in the
  # middle of the omega M-step -- after the covstruct mask but BEFORE the
  # variance floors, the fix()ed-value restore and the diagonal enforcement.
  # The kernel sampled with the fixed value while `fit$omega` reported the
  # unconstrained M-step estimate.
  .mod <- function() {
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

  .ctl <- saemControl(nBurn = 20, nEm = 20, print = 0, seed = 42,
                      calcTables = FALSE, covMethod = "")

  .fit <- suppressMessages(suppressWarnings(
    nlmixr2(.mod, nlmixr2data::theo_sd, "saem", .ctl)))

  test_that("saem reports a fix()ed eta variance as the fixed value (#1073)", {
    expect_equal(.fit$omega["eta.ka", "eta.ka"], 0.3, tolerance = 1e-12)
    # the unfixed variance is still estimated away from its ini()
    expect_true(.fit$omega["eta.cl", "eta.cl"] != 0.1)
  })

  test_that("the reporting snapshot is taken after the constraints (#1073)", {
    # mechanism: with a single component nothing is pooled, so the reported
    # matrix must be the live one the sampler used -- that equality is exactly
    # what the mis-ordered snapshot broke.
    expect_equal(.fit$saem$Gamma2_phi1Report, .fit$saem$Gamma2_phi1,
                 tolerance = 1e-12)
    expect_equal(.fit$saem$Gamma2_phi1Report[1, 1], 0.3, tolerance = 1e-12)
  })

  # A non-mu-referenced eta owns a phi column the same way, so it drifted the
  # same way.
  .modNonMu <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ fix(1)
      eta.cl ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka) * exp(eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("a non-mu-referenced fix()ed eta variance reports the fixed value (#1073)", {
    .fitNonMu <- suppressMessages(suppressWarnings(
      nlmixr2(.modNonMu, nlmixr2data::theo_sd, "saem", .ctl)))
    expect_equal(.fitNonMu$omega["eta.ka", "eta.ka"], 1, tolerance = 1e-12)
  })
})
